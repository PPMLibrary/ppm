      minclude ppm_header(ppm_module_name="loadbal_map_subpart")

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE loadbal_map_subpart_s(topoid,Pc,isSender,isub,from,to,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE loadbal_map_subpart_d(topoid,Pc,isSender,isub,from,to,info)
#endif
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
      USE ppm_module_map
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
      !!! ID of the topology
#if    __KIND == __SINGLE_PRECISION
      TYPE(ppm_t_particles_s) , INTENT(INOUT) :: Pc
#elif __KIND == __DOUBLE_PRECISION
      TYPE(ppm_t_particles_d) , INTENT(INOUT) :: Pc
#endif
      LOGICAL                 , INTENT(IN   ) :: isSender
      !!! Am I the sender of isub?
      INTEGER                 , INTENT(IN   ) :: isub
      !!! Global ID of the subdomain to be sent

      INTEGER                 , INTENT(IN   ) :: from
      !!! MPI rank of the sender
      INTEGER                 , INTENT(IN   ) :: to
      !!! MPI rank of the receiver
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)          :: ldu
      INTEGER, DIMENSION(:), POINTER :: bcdef => NULL(),idlbcommseq=>NULL()
      INTEGER                        :: ndlbcommseq
      INTEGER                        :: i,j,k,idom,ipart,nlist1,nlist2
      INTEGER                        :: sendrank,recvrank,ltopoid
      INTEGER                        :: nneighsubs, jdom, local_isub
      INTEGER                        :: iopt,iset,ibuffer,isonneigh
      INTEGER                        :: recvidx,Npart,Npart_new,maxneigh
      CHARACTER(ppm_char)            :: mesg
      REAL(MK)                       :: t0
      LOGICAL                        :: valid
      LOGICAL                        :: ignoreunassigned
      TYPE(ppm_t_topo)    , POINTER  :: topo => NULL()
      REAL(MK), DIMENSION(:,:), POINTER     :: xp

      !----------------------------------------------------------------------
      !  Work lists
      !----------------------------------------------------------------------
!      INTEGER, DIMENSION(:), POINTER :: ilist1     => NULL()
!      INTEGER, DIMENSION(:), POINTER :: ilist2     => NULL()
!      INTEGER, DIMENSION(:), POINTER :: part2proc  => NULL()
!      INTEGER, DIMENSION(:), POINTER :: ineighsubs => NULL()
#if    __KIND == __SINGLE_PRECISION
      CLASS(ppm_t_part_prop_s_), POINTER :: prop => NULL()
      CLASS(ppm_t_neighlist_s_), POINTER :: nl => NULL()
      CLASS(ppm_t_operator_discr_),   POINTER :: op => NULL()
#elif __KIND == __DOUBLE_PRECISION
      CLASS(ppm_t_part_prop_d_), POINTER :: prop => NULL()
      CLASS(ppm_t_neighlist_d_), POINTER :: nl => NULL()
      CLASS(ppm_t_operator_discr_),   POINTER :: op => NULL()
#endif
      !   PRIVATE :: ilist1,ilist2,part2proc,ineighsubs
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      start_subroutine("loadbal_map_subpart")
      topo => ppm_topo(topoid)%t
      !-------------------------------------------------------------------------
      !  isub will be sent to the underloaded processor
      !-------------------------------------------------------------------------

      bcdef => topo%bcdef
      Npart = Pc%Npart
!      stdout("...Npart",Npart)
      xp => Pc%xp
      ltopoid = Pc%active_topoid
      ndlbcommseq = 2
!      stdout("from-beginning",from,"to",to)
      !-------------------------------------------------------------------------
      !  Impose the boundary conditions.
      !  This is now taken care of transparently. Unless the user needs to
      !  impose the boundary conditions at a very specific stage of his
      !  simulation there is no need to bother the user with this.
      !-------------------------------------------------------------------------
      CALL ppm_impose_part_bc(topoid,xp,Npart,info)
      or_fail("imposing particle BCs failed.")
!      stdout("...Npart after impose",Npart)
      ! if there is still some data left in the buffer, warn the user
      IF (ppm_buffer_set .GT. 0) THEN
        fail("Buffer was not empty. Possible loss of data!")
      ENDIF

      ! first check if the optimal communication protocol is known
      IF (.NOT. topo%isoptimized) THEN
        ! if not: determine it before calling map_part_partial
        CALL ppm_util_commopt(topoid,info)
        IF (info.NE.0) GOTO 9999
        IF (ppm_debug .GT. 1) THEN
            DO i=1,topo%nneighproc
                WRITE(mesg,'(A,I4)') 'have neighbor: ',  topo%ineighproc(i)
                CALL ppm_write(ppm_rank,'ppm_loadbal_map_subpart',mesg,info)
            END DO
            DO i=1,topo%ncommseq
                WRITE(mesg,'(A,I4)') 'communicate: ', topo%icommseq(i)
                CALL ppm_write(ppm_rank,'ppm_loadbal_map_subpart',mesg,info)
            END DO
        ENDIF
      END IF

      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_partial

      !-------------------------------------------------------------------------
      !  Allocate memory for particle lists & comm. sequence for dlb
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = Npart
      CALL ppm_alloc(ilist1,ldu,iopt,info)
       or_fail_alloc("ilist1")
      CALL ppm_alloc(ilist2,ldu,iopt,info)
       or_fail_alloc("ilist2")
      ldu(1) = ndlbcommseq
      CALL ppm_alloc(idlbcommseq,ldu,iopt,info)
       or_fail_alloc("idlbcommseq")
!      stdout("topo_commseq",'topo%icommseq')
      !-------------------------------------------------------------------------
      !  initialize the idlbcommseq array
      !-------------------------------------------------------------------------
      idlbcommseq(1) = ppm_rank
      IF (isSender) THEN
        idlbcommseq(2) = to
      ELSE
        idlbcommseq(2) = from
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for neighbor list (subs)
      !  This could be optimized by dynamically growing the array whenever
      !  a neighbor is found. Right now, I just alloc the maximum number
      !  possible (nsubs is not too big normally).
      !-------------------------------------------------------------------------
      ldu(1) = topo%nsubs
      CALL ppm_alloc(ineighsubs,ldu,iopt,info)
       or_fail_alloc("neighbor subs list INEIGHSUBS")

      !-------------------------------------------------------------------------
      !  Find all subs residing on neighboring processors (or on local one)
      !  NEED TO INCLUDE THE LOCAL SUBS AS WELL because otherwise:
      !    (1) we will not be able to assign all particles and the check will
      !        fail
      !    (2) The ppm_map_part_send routine expects the first block in the
      !        buffer to be the local particles and SKIPS IT.
      !-------------------------------------------------------------------------
      nneighsubs = 0
!      stdout("from",from,"to",to)
      topo%sub2proc(isub)  = to

      DO i=1,topo%nsubs
          ! if either on local or the receiving proc
          IF (topo%sub2proc(i) .EQ. from .OR.topo%sub2proc(i).EQ.to) THEN
              nneighsubs             = nneighsubs + 1
              ineighsubs(nneighsubs) = i
          ENDIF
      ENDDO

!      print*,ppm_rank,'ineighsubs',ineighsubs
!      print*,ppm_rank,'nneighsubs',nneighsubs
      !-------------------------------------------------------------------------
      !  Allocate memory for the pointer to the buffer; for the partial map we
      !  need entries for each communication round. Thus ldu(1) =
      !  ncommseq(topoid) + 1
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ndlbcommseq + 1
      CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
       or_fail_alloc("PPM_PSENDBUFFER")
      !-------------------------------------------------------------------------
      !  Alloc particle to processor and buffer to particle maps
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = Npart
      CALL ppm_alloc(part2proc,ldu,iopt,info)
       or_fail_alloc("part2proc")
      ldu(1) = Npart
      CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
       or_fail_alloc("ppm_buffer2part")

      !-------------------------------------------------------------------------
      !  Initialize the particle list and part2proc list
      !-------------------------------------------------------------------------
!      stdout("...Npart..before nlist1",Npart)
      nlist1 = Npart
      DO ipart=1,Npart
         ilist1(ipart)    = ipart
         part2proc(ipart) = -1
      ENDDO

      !-------------------------------------------------------------------------
      !  Assign particles to processors in part2proc(ipart) = sendrank
      !-------------------------------------------------------------------------
      IF (ppm_dim.EQ.2) THEN
         !----------------------------------------------------------------------
         !  Loop over the subdomains on the neighboring processors (2D)
         !----------------------------------------------------------------------
         DO jdom=nneighsubs,1,-1
            idom = ineighsubs(jdom)
            sendrank = topo%sub2proc(idom)

            !-------------------------------------------------------------------
            !  Loop over the remaining particles not yet assigned to a processor
            !-------------------------------------------------------------------
            nlist2 = 0
            DO i=1,nlist1
               ipart = ilist1(i)
               !----------------------------------------------------------------
               !  If the particle is inside the current subdomain, assign it
               !----------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
               IF (xp(1,ipart) .GE. topo%min_subs(1,idom).AND.   &
     &             xp(1,ipart) .LE. topo%max_subs(1,idom).AND.   &
     &             xp(2,ipart) .GE. topo%min_subs(2,idom).AND.   &
     &             xp(2,ipart) .LE. topo%max_subs(2,idom)) THEN
                   !------------------------------------------------------------
                   !  In the non-periodic case, allow particles that are
                   !  exactly ON an upper EXTERNAL boundary.
                   !------------------------------------------------------------
                   IF((xp(1,ipart).LT.topo%max_subs(1,idom) .OR.  &
     &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
     &                (xp(2,ipart).LT.topo%max_subs(2,idom) .OR.  &
     &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
     &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
#elif  __KIND == __DOUBLE_PRECISION
               IF (xp(1,ipart) .GE. topo%min_subd(1,idom).AND.   &
     &             xp(1,ipart) .LE. topo%max_subd(1,idom).AND.   &
     &             xp(2,ipart) .GE. topo%min_subd(2,idom).AND.   &
     &             xp(2,ipart) .LE. topo%max_subd(2,idom)) THEN
                   !------------------------------------------------------------
                   !  In the non-periodic case, allow particles that are
                   !  exactly ON an upper EXTERNAL boundary.
                   !------------------------------------------------------------
                   IF((xp(1,ipart).LT.topo%max_subd(1,idom) .OR.  &
     &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
     &                (xp(2,ipart).LT.topo%max_subd(2,idom) .OR.  &
     &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
     &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
#endif
                       part2proc(ipart) = sendrank
                   ELSE
                      !---------------------------------------------------------
                      !  if not, add it the list of particle to search
                      !---------------------------------------------------------
                      nlist2         = nlist2 + 1
                      ilist2(nlist2) = ipart
                   ENDIF
                ELSE
                   !---------------------------------------------------------
                   !  if not, add it the list of particle to search
                   !---------------------------------------------------------
                   nlist2         = nlist2 + 1
                   ilist2(nlist2) = ipart
                ENDIF
            ENDDO

            !-------------------------------------------------------------------
            !  Copy the lists (well, only if nlist2 changed - decreased)
            !-------------------------------------------------------------------
            IF (nlist2.NE.nlist1) THEN
               nlist1 = nlist2
               DO i=1,nlist1
                  ilist1(i) = ilist2(i)
               ENDDO
            ENDIF

            !-------------------------------------------------------------------
            !  Exit if the list is empty
            !-------------------------------------------------------------------
            IF (nlist1.EQ.0) EXIT
         ENDDO
      ELSE
         !----------------------------------------------------------------------
         !  Loop over the subdomains (since the first domains are most likely
         !  to be empty, we look backwards to reduce the number of elements in
         !  nlist2 as fast as possible) (3D)
         !----------------------------------------------------------------------
         DO jdom=nneighsubs,1,-1
            idom = ineighsubs(jdom)
            sendrank   = topo%sub2proc(idom)
            !-------------------------------------------------------------------
            !  Loop over the remaining particles
            !-------------------------------------------------------------------
            nlist2 = 0
            DO i=1,nlist1
               ipart = ilist1(i)
               !----------------------------------------------------------------
               !  If the particle is inside the current subdomain, assign it
               !----------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
               IF (xp(1,ipart) .GE. topo%min_subs(1,idom).AND.   &
     &             xp(1,ipart) .LE. topo%max_subs(1,idom).AND.   &
     &             xp(2,ipart) .GE. topo%min_subs(2,idom).AND.   &
     &             xp(2,ipart) .LE. topo%max_subs(2,idom).AND.   &
     &             xp(3,ipart) .GE. topo%min_subs(3,idom).AND.   &
     &             xp(3,ipart) .LE. topo%max_subs(3,idom)) THEN
                   !------------------------------------------------------------
                   !  In the non-periodic case, allow particles that are
                   !  exactly ON an upper EXTERNAL boundary.
                   !------------------------------------------------------------
                   IF((xp(1,ipart) .LT. topo%max_subs(1,idom) .OR. &
     &                (topo%subs_bc(2,idom).EQ.1           .AND. &
     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND. &
     &                (xp(2,ipart) .LT. topo%max_subs(2,idom) .OR. &
     &                (topo%subs_bc(4,idom).EQ.1           .AND. &
     &                bcdef(4) .NE. ppm_param_bcdef_periodic))    .AND. &
     &                (xp(3,ipart) .LT. topo%max_subs(3,idom) .OR. &
     &                (topo%subs_bc(6,idom).EQ.1           .AND. &
     &                bcdef(6) .NE. ppm_param_bcdef_periodic))   ) THEN
#elif  __KIND == __DOUBLE_PRECISION
               IF (xp(1,ipart) .GE. topo%min_subd(1,idom).AND.   &
     &             xp(1,ipart) .LE. topo%max_subd(1,idom).AND.   &
     &             xp(2,ipart) .GE. topo%min_subd(2,idom).AND.   &
     &             xp(2,ipart) .LE. topo%max_subd(2,idom).AND.   &
     &             xp(3,ipart) .GE. topo%min_subd(3,idom).AND.   &
     &             xp(3,ipart) .LE. topo%max_subd(3,idom)) THEN
                   !------------------------------------------------------------
                   !  In the non-periodic case, allow particles that are
                   !  exactly ON an upper EXTERNAL boundary.
                   !------------------------------------------------------------
                   IF((xp(1,ipart).LT.topo%max_subd(1,idom) .OR. &
     &                (topo%subs_bc(2,idom).EQ.1           .AND. &
     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND. &
     &                (xp(2,ipart).LT.topo%max_subd(2,idom) .OR. &
     &                (topo%subs_bc(4,idom).EQ.1           .AND. &
     &                bcdef(4).NE. ppm_param_bcdef_periodic))    .AND. &
     &                (xp(3,ipart).LT.topo%max_subd(3,idom) .OR. &
     &                (topo%subs_bc(6,idom).EQ.1           .AND. &
     &                bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN
#endif
                      part2proc(ipart) = sendrank
                   ELSE
                      nlist2         = nlist2 + 1
                      ilist2(nlist2) = ipart
                   ENDIF
                ELSE
                   nlist2         = nlist2 + 1
                   ilist2(nlist2) = ipart
                ENDIF
            ENDDO

            !-------------------------------------------------------------------
            !  Copy the lists (well, only if nlist2 changed - decreased)
            !-------------------------------------------------------------------
            IF (nlist2.NE.nlist1) THEN
               nlist1 = nlist2
               DO i=1,nlist1
                  ilist1(i) = ilist2(i)
               ENDDO
            ENDIF

            !-------------------------------------------------------------------
            !  Exit if the list is empty
            !-------------------------------------------------------------------
            IF (nlist1.EQ.0) EXIT
         ENDDO
      ENDIF
!
      !-------------------------------------------------------------------------
      !  Check if we sold all the particles. If not some of them have move
      !  too far we return an info and the user should call the global map.
      !-------------------------------------------------------------------------
!      IF (.NOT. ignoreunassigned) THEN
!         IF (nlist2.GT.0) THEN
!            stdout("ERROR::nlist2>0")
!            fail("Please call ppm_map_part_global")
!         ENDIF
!      ENDIF

      !-------------------------------------------------------------------------
      !  Store the number of buffer entries (this is the first)
      !-------------------------------------------------------------------------
      ppm_buffer_set = 1

      !-------------------------------------------------------------------------
      !  Allocate memory for the send buffer meta data
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_buffer_set
      CALL ppm_alloc(ppm_buffer_dim ,ldu,iopt,info)
       or_fail_alloc("ppm_buffer_dim")
      CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
       or_fail_alloc("ppm_buffer_type")

      ppm_buffer_dim(ppm_buffer_set)  = ppm_dim
#if    __KIND == __SINGLE_PRECISION
      ppm_buffer_type(ppm_buffer_set) = ppm_kind_single
#else
      ppm_buffer_type(ppm_buffer_set) = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  (Re)allocate memory for the buffer
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_dim*Npart
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
      ELSE
         CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
      ENDIF
       or_fail_alloc("PPM_SENDBUFFER")

      !-------------------------------------------------------------------------
      !  Allocate memory for the sendlist
      !-------------------------------------------------------------------------
      ppm_nsendlist = ndlbcommseq
      ppm_nrecvlist = ndlbcommseq
      ldu(1)        = ppm_nsendlist
      CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
       or_fail_alloc("ppm_isendlist")
      CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
       or_fail_alloc("ppm_irecvlist")

      !-------------------------------------------------------------------------
      !  Initialize the particle lists
      !-------------------------------------------------------------------------
!      stdout("...Npart...after nlist1 (2)",Npart)
      nlist1 = Npart
      DO ipart=1,Npart
         ilist1(ipart) = ipart
      ENDDO

      !-------------------------------------------------------------------------
      !  loop over the neighboring processors according to the optimized
      !  communication sequence.
      !-------------------------------------------------------------------------
      ppm_psendbuffer(1) = 1
      ppm_nsendlist      = 0
      ppm_nrecvlist      = 0
      iset               = 0
      ibuffer            = 0

!      ppm_isendlist(ppm_nsendlist) = to
!      ppm_irecvlist(ppm_nrecvlist) = from
      DO i=1,ndlbcommseq
         !----------------------------------------------------------------------
         !  get next neighbor to send/recv to/from.
         !----------------------------------------------------------------------
         sendrank = idlbcommseq(i)
         recvrank = sendrank

         !----------------------------------------------------------------------
         !  Store the processor to which we will send to
         !----------------------------------------------------------------------
         ppm_nsendlist                = ppm_nsendlist + 1
         ppm_isendlist(ppm_nsendlist) = sendrank

         !----------------------------------------------------------------------
         !  Store the processor to which we will recv from
         !----------------------------------------------------------------------
         ppm_nrecvlist                = ppm_nrecvlist + 1
         ppm_irecvlist(ppm_nrecvlist) = recvrank

         !----------------------------------------------------------------------
         !  Only assign particles if there is any communication for this
         !  processor in this round
         !----------------------------------------------------------------------
         IF (sendrank .GE. 0) THEN

             !------------------------------------------------------------------
             !  Initialize the buffer count
             !------------------------------------------------------------------
             nlist2 = 0
             IF (ppm_dim .EQ. 2) THEN
                DO j=1,nlist1
                   ipart = ilist1(j)
                   IF (part2proc(ipart).EQ.sendrank) THEN
                      !---------------------------------------------------------
                      !  increment the buffer counter
                      !---------------------------------------------------------
                      iset = iset + 1

                      !---------------------------------------------------------
                      !  Store the id of the particle
                      !---------------------------------------------------------
                      ppm_buffer2part(iset) = ipart

                      !---------------------------------------------------------
                      !  Store the particle
                      !---------------------------------------------------------
                      IF (ppm_kind .EQ. ppm_kind_double) THEN
#if    __KIND == __SINGLE_PRECISION
                         ibuffer                  = ibuffer + 1
                         ppm_sendbufferd(ibuffer) = REAL(xp(1,ipart),   &
     &                       ppm_kind_double)
                         ibuffer                  = ibuffer + 1
                         ppm_sendbufferd(ibuffer) = REAL(xp(2,ipart),   &
     &                       ppm_kind_double)
#else
                         ibuffer                  = ibuffer + 1
                         ppm_sendbufferd(ibuffer) = xp(1,ipart)
                         ibuffer                  = ibuffer + 1
                         ppm_sendbufferd(ibuffer) = xp(2,ipart)
#endif
                      ELSE
#if    __KIND == __SINGLE_PRECISION
                         ibuffer                  = ibuffer + 1
                         ppm_sendbuffers(ibuffer) = xp(1,ipart)
                         ibuffer                  = ibuffer + 1
                         ppm_sendbuffers(ibuffer) = xp(2,ipart)
#else
                         ibuffer                  = ibuffer + 1
                         ppm_sendbuffers(ibuffer) = REAL(xp(1,ipart),   &
     &                       ppm_kind_single)
                         ibuffer                  = ibuffer + 1
                         ppm_sendbuffers(ibuffer) = REAL(xp(2,ipart),   &
     &                       ppm_kind_single)
#endif
                      ENDIF
                   ELSE
                      nlist2         = nlist2 + 1
                      ilist2(nlist2) = ipart
                   ENDIF
                ENDDO
             ELSE
                DO j=1,nlist1
                   ipart = ilist1(j)
                   IF (part2proc(ipart).EQ.sendrank) THEN
                      !---------------------------------------------------------
                      !  increment the buffer counter
                      !---------------------------------------------------------
                      iset = iset + 1

                      !---------------------------------------------------------
                      !  Store the id of the particle
                      !---------------------------------------------------------
                      ppm_buffer2part(iset) = ipart

                      !---------------------------------------------------------
                      !  Store the particle
                      !---------------------------------------------------------
                      IF (ppm_kind .EQ. ppm_kind_double) THEN
#if    __KIND == __SINGLE_PRECISION
                         ibuffer                  = ibuffer + 1
                         ppm_sendbufferd(ibuffer) = REAL(xp(1,ipart),   &
     &                       ppm_kind_double)
                         ibuffer                  = ibuffer + 1
                         ppm_sendbufferd(ibuffer) = REAL(xp(2,ipart),   &
     &                       ppm_kind_double)
                         ibuffer                  = ibuffer + 1
                         ppm_sendbufferd(ibuffer) = REAL(xp(3,ipart),   &
     &                       ppm_kind_double)
#else
                         ibuffer                  = ibuffer + 1
                         ppm_sendbufferd(ibuffer) = xp(1,ipart)
                         ibuffer                  = ibuffer + 1
                         ppm_sendbufferd(ibuffer) = xp(2,ipart)
                         ibuffer                  = ibuffer + 1
                         ppm_sendbufferd(ibuffer) = xp(3,ipart)
#endif
                      ELSE
#if    __KIND == __SINGLE_PRECISION
                         ibuffer                  = ibuffer + 1
                         ppm_sendbuffers(ibuffer) = xp(1,ipart)
                         ibuffer                  = ibuffer + 1
                         ppm_sendbuffers(ibuffer) = xp(2,ipart)
                         ibuffer                  = ibuffer + 1
                         ppm_sendbuffers(ibuffer) = xp(3,ipart)
#else
                         ibuffer                  = ibuffer + 1
                         ppm_sendbuffers(ibuffer) = REAL(xp(1,ipart),   &
     &                       ppm_kind_single)
                         ibuffer                  = ibuffer + 1
                         ppm_sendbuffers(ibuffer) = REAL(xp(2,ipart),   &
     &                       ppm_kind_single)
                         ibuffer                  = ibuffer + 1
                         ppm_sendbuffers(ibuffer) = REAL(xp(3,ipart),   &
     &                       ppm_kind_single)
#endif
                      ENDIF
                   ELSE
                      nlist2         = nlist2 + 1
                      ilist2(nlist2) = ipart
                   ENDIF
                ENDDO
             ENDIF

             !------------------------------------------------------------------
             !  Swap the lists
             !------------------------------------------------------------------
             nlist1 = nlist2
             DO j=1,nlist1
                ilist1(j) = ilist2(j)
             ENDDO
         ENDIF     ! sendrank .GE. 0

         !----------------------------------------------------------------------
         !  Update the buffer pointer
         !----------------------------------------------------------------------
         ppm_psendbuffer(i+1) = iset + 1

      ENDDO

      !-------------------------------------------------------------------------
      !  All particles have to go
      !-------------------------------------------------------------------------
      ppm_psendbuffer(ndlbcommseq+1) = Npart + 1 - nlist2
!      ppm_psendbuffer(ndlbcommseq+1) = ppm_psendbuffer(ndlbcommseq)
      !-------------------------------------------------------------------------
      !  Store the current size of the buffer
      !-------------------------------------------------------------------------
      ppm_nsendbuffer = ibuffer

!      stdout("LOADBAL:: ppm_psendbuffer:",'ppm_psendbuffer(:)',"ppm_nsendbuffer",ppm_nsendbuffer)
      !-------------------------------------------------------------------------
      !  Deallocate the memory for the lists
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(ilist1,ldu,iopt,info)
      or_fail_alloc("ilist1")
      CALL ppm_alloc(ilist2,ldu,iopt,info)
      or_fail_alloc("ilist2")
      CALL ppm_alloc(part2proc,ldu,iopt,info)
      or_fail_alloc("part2proc")
      CALL ppm_alloc(ineighsubs,ldu,iopt,info)
      or_fail_alloc("ineighsubs")
      Pc%stats%nb_part_map = Pc%stats%nb_part_map + 1
!      nlist2 = 0
!      stdout("Pc-Npart",'Pc%Npart',"Npart",Npart)
!      Pc%Npart = Npart
      !-------------------------------------------------------------------------
      !  Buffers are filled and ready to be sent
      !-------------------------------------------------------------------------
!      stdout("buffers are full!")
      prop => Pc%props%begin()
      DO WHILE (ASSOCIATED(prop))
          IF (prop%flags(ppm_ppt_map_parts)) THEN
              CALL Pc%map_part_push_legacy(Pc%props%iter_id,info)
              or_fail("ppm_map_part_push_legacy")
          ENDIF
          prop => Pc%props%next()
      ENDDO
!      stdout("Npart before map_part_send",'Pc%Npart')
      CALL ppm_map_part_send(Pc%Npart,Npart_new,info)
      or_fail("ppm_map_part_send")
!      stdout("...Npart AND _new",'Pc%Npart',Npart_new)
!      stdout("LOADBAL2:: ppm_nsendlist:",ppm_nsendlist,"ppm_nsendbuffer",ppm_nsendbuffer)
      prop => Pc%props%last()
      DO WHILE (ASSOCIATED(prop))
          IF (prop%flags(ppm_ppt_map_parts)) THEN
              CALL Pc%map_part_pop_legacy(Pc%props%iter_id,Npart_new,info)
              or_fail("ppm_map_part_pop_legacy")
              prop%flags(ppm_ppt_partial) = .TRUE.
          ENDIF
          prop => Pc%props%prev()
      ENDDO

      CALL ppm_map_part_pop(Pc%xp,ppm_dim,Pc%Npart,Npart_new,info)
      or_fail("ppm_map_part_prop")

      ! Update states
      ! Number of particles on this processor
!      stdout("...Npart_new",Npart_new)
      Pc%Npart = Npart_new
      Pc%Mpart = Pc%Npart

      ! Particles are now mapped on the active topology
      Pc%flags(ppm_part_partial) = .TRUE.
      ! Pc have been re-indexed and ghosts have not been computed
      Pc%flags(ppm_part_ghosts) = .FALSE.

      ! values for poperty arrays have been mapped and ghosts
      ! are no longer up-to-date
      prop => Pc%props%begin()
      DO WHILE (ASSOCIATED(prop))
          prop%flags(ppm_ppt_ghosts) = .FALSE.
          prop => Pc%props%next()
      ENDDO

      ! particles have been re-indexed and neighbour lists not updated
      nl => Pc%neighs%begin()
      DO WHILE (ASSOCIATED(nl))
          nl%uptodate = .FALSE.
          nl => Pc%neighs%next()
      ENDDO
      Pc%flags(ppm_part_neighlists) = .FALSE.

      ! particles have been re-indexed and operators need be recomputed
      IF (ASSOCIATED(Pc%ops)) THEN
          op => Pc%ops%begin()
          DO WHILE (ASSOCIATED(op))
              op%flags(ppm_ops_iscomputed) = .FALSE.
              op => Pc%ops%next()
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Update topology related info
      !-------------------------------------------------------------------------
!      print*,ppm_rank,'isublist BEFORE:',topo%isublist
      IF (isSender) THEN

!        stdout("global subID to be sent",isub)
        DO i=1,topo%nsublist
            IF (topo%isublist(i).EQ.isub) THEN
!                stdout("local subID is found")
                local_isub = i
            ENDIF
        ENDDO
        maxneigh = SIZE(topo%ineighsubs(:,local_isub))
        stdout("maxneigh",maxneigh)
        CALL ppm_loadbal_sendsub(topoid,local_isub,to,maxneigh,topo%prec, &
     &                           info)
        or_fail("ppm_loadbal_sendsub")
      ELSE
!        stdout("from end",from,"to",to)
        CALL ppm_loadbal_recvsub(topoid,isub,from,topo%prec,info)
        or_fail("ppm_loadbal_recvsub")
      ENDIF
!      stdout("send/recv complete!!")
      ! update the ghosts too
!      CALL Pc%map_ghosts(info)
!      or_fail("map_ghosts in loadbal failed")

!      print*,ppm_rank,'sub2proc:',topo%sub2proc
!      stdout("sub2proc",'topo%sub2proc')

!      print*,ppm_rank,'isublist AFTER:',topo%isublist
!      stdout("isublist",'topo%isublist')


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
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE loadbal_map_subpart_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE loadbal_map_subpart_d
#endif


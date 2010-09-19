      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_map_part_partial
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Lab (ETH Zurich), 
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

#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_partial_s(topoid,xp,Npart,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_partial_d(topoid,xp,Npart,info)
      !!! This routine maps the particles onto the topology using a local map
      !!! (i.e. each processor only communicates with its neighbors).
      !!!
      !!! [NOTE]
      !!! The first part of the buffer contains the on processor data.

#endif 
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_check_id
      USE ppm_module_typedef
      USE ppm_module_write
      USE ppm_module_util_commopt
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! Particle coordinates
      INTEGER                 , INTENT(IN   ) :: Npart
      !!! The number of particles (on the local processor)
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! ID of the topology
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------

      INTEGER, DIMENSION(3)          :: ldu
      INTEGER, DIMENSION(:), POINTER :: bcdef
      INTEGER                        :: i,j,k,idom,ipart,nlist1,nlist2
      INTEGER                        :: sendrank,recvrank
      INTEGER                        :: nneighsubs, jdom
      INTEGER                        :: iopt,iset,ibuffer,isonneigh
      INTEGER                        :: recvidx
      CHARACTER(ppm_char)            :: mesg
      REAL(MK)                       :: t0
      LOGICAL                        :: valid
      TYPE(ppm_t_topo)    , POINTER  :: topo
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_partial',t0,info)



      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t
      bcdef => topo%bcdef


      ! if there is still some data left in the buffer, warn the user
      IF (ppm_buffer_set .GT. 0) THEN
        info = ppm_error_warning
        CALL ppm_error(ppm_err_map_incomp,'ppm_map_part_partial',  &
     &      'Buffer was not empty. Possible loss of data!',__LINE__,info)
      ENDIF


      ! first check if the optimal communication protocol is known
      IF (.NOT. topo%isoptimized) THEN
        ! if not: determine it before calling map_part_partial
        CALL ppm_util_commopt(topoid,info)
        IF (info.NE.0) GOTO 9999
        IF (ppm_debug .GT. 1) THEN
            DO i=1,topo%nneighproc
                WRITE(mesg,'(A,I4)') 'have neighbor: ',  topo%ineighproc(i)
                CALL ppm_write(ppm_rank,'ppm_map_part_partial',mesg,info)
            END DO
            DO i=1,topo%ncommseq
                WRITE(mesg,'(A,I4)') 'communicate: ', topo%icommseq(i)
                CALL ppm_write(ppm_rank,'ppm_map_part_partial',mesg,info)
            END DO
        ENDIF
      END IF

      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls 
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_partial

      !-------------------------------------------------------------------------
      !  Allocate memory for particle lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = Npart
      CALL ppm_alloc(ilist1,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_partial',     &
     &        'particle list 1 ILIST1',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_partial',     &
     &        'particle list 2 ILIST2',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for neighbor list (subs)
      !  This could be optimized by dynamically growing the array whenever
      !  a neighbor is found. Right now, I just alloc the maximum number
      !  possible (nsubs is not too big normally).
      !-------------------------------------------------------------------------
      ldu(1) = topo%nsubs
      CALL ppm_alloc(ineighsubs,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_partial',     &
     &        'neighbor subs list INEIGHSUBS',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find all subs residing on neighboring processors (or on local one)
      !  NEED TO INCLUDE THE LOCAL SUBS AS WELL because otherwise:
      !    (1) we will not be able to assign all particles and the check will
      !        fail
      !    (2) The ppm_map_part_send routine expects the first block in the
      !        buffer to be the local particles and SKIPS IT.
      !-------------------------------------------------------------------------
      nneighsubs = 0
      DO i=1,topo%nsubs
          ! if either on local ...
          IF (topo%sub2proc(i) .EQ. ppm_rank) THEN
              nneighsubs             = nneighsubs + 1
              ineighsubs(nneighsubs) = i
          ELSE
              ! or on any of my neighbors
              isonneigh = 0
              j = 1
              DO WHILE (isonneigh .EQ. 0 .AND. j .LE. topo%nneighproc)
                  IF (topo%sub2proc(i) .EQ. topo%ineighproc(j)) THEN
                      nneighsubs             = nneighsubs + 1
                      ineighsubs(nneighsubs) = i
                      isonneigh = 1
                  ENDIF
                  j = j + 1
              ENDDO
          ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for the pointer to the buffer; for the partial map we
      !  need entries for each communication round. Thus ldu(1) = 
      !  ncommseq(topoid) + 1
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = topo%ncommseq + 1
      CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_partial',     &
     &        'particle send buffer PPM_PSENDBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Alloc particle to processor and buffer to particle maps
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = Npart
      CALL ppm_alloc(part2proc,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_partial',     &
     &        'particles-to-processors map PART2PROC',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = Npart
      CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_partial',     &
     &        'buffer-to-particles map PPM_BUFFER2PART',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the particle list and part2proc list
      !-------------------------------------------------------------------------
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

      !-------------------------------------------------------------------------
      !  Check if we sold all the particles. If not some of them have move
      !  too far we return an info and the user should call the global map.
      !-------------------------------------------------------------------------
      IF (nlist2.GT.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_part_unass,'ppm_map_part_partial', &
     &       'Please call ppm_map_part_global',__LINE__,info)
      ENDIF 

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
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_partial',     &
     &        'buffer dimensions PPM_BUFFER_DIM',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_partial',     &
     &        'buffer types PPM_BUFFER_TYPE',__LINE__,info)
          GOTO 9999
      ENDIF

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
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_partial',     &
     &        'global send buffer PPM_SENDBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the sendlist
      !-------------------------------------------------------------------------
      ppm_nsendlist = topo%ncommseq
      ppm_nrecvlist = topo%ncommseq
      ldu(1)        = ppm_nsendlist
      CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_partial',     &
     &        'send list PPM_ISENDLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_partial',     &
     &        'receive list PPM_IRECVLIST',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the particle lists 
      !-------------------------------------------------------------------------
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
      DO i=1,topo%ncommseq
         !----------------------------------------------------------------------
         !  get next neighbor to send/recv to/from.
         !----------------------------------------------------------------------
         sendrank = topo%icommseq(i)
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
      ppm_psendbuffer(topo%ncommseq+1) = npart + 1

      !-------------------------------------------------------------------------
      !  Store the current size of the buffer
      !-------------------------------------------------------------------------
      ppm_nsendbuffer = ibuffer

      !-------------------------------------------------------------------------
      !  Deallocate the memory for the lists
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(ilist1,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_partial',     &
     &        'particle list 1 ILIST1',__LINE__,info)
      ENDIF
      CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_partial',     &
     &        'particle list 2 ILIST2',__LINE__,info)
      ENDIF
      CALL ppm_alloc(part2proc,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_partial',     &
     &        'particles-to-processors map PART2PROC',__LINE__,info)
      ENDIF
      CALL ppm_alloc(ineighsubs,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_partial',     &
     &        'neighbor subs list INEIGHSUBS',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_partial',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_partial',  &
     &            'Npart must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (topoid .EQ. ppm_param_topo_undefined) THEN
            info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_partial',  &
     &            'topoid must be the ID of a defined topology',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_partial',  &
     &             'topoid out of range',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_partial_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_partial_d
#endif

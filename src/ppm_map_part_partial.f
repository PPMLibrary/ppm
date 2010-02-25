      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_map_part_partial
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine maps the particles onto the current 
      !                 topology using a local map (i.e. each processor
      !                 only communicates with its neighbors).
      !
      !  Input        : xp(:,:)    (F) the position of the particles
      !                 Npart      (I) the number of particles (on the
      !                                local processor)
      !
      !  Input/output :
      !
      !  Output       : info       (I) return status. 0 on success.
      !
      !  Remarks      : The first part of the buffer contains the on processor
      !                 data. No topoid as argument: always use the current
      !                 topology as this is just an update and using it on
      !                 a topology change would badly mess up tings.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_part_partial.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.24  2006/09/04 18:34:51  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.23  2006/02/06 15:29:10  michaebe
      !  fixed bug 0000016
      !
      !  Revision 1.22  2005/11/30 12:59:50  ivos
      !  Boundary conditions are now explicitly accounted for. Particles
      !  exactly ON the upper domain boundary are allowed in all
      !  non-periodic cases.
      !
      !  Revision 1.21  2004/11/11 15:26:17  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.20  2004/10/01 16:09:07  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.19  2004/08/30 07:49:18  ivos
      !  Relaxed argument check from Npart.LE.0 to Npart.LT.0.
      !
      !  Revision 1.18  2004/08/03 12:48:28  ivos
      !  Explicitly split loops for ppm_dim.EQ.2 and ppm_dim.EQ.3 cases.
      !  All loops now vectorize on the NEC.
      !
      !  Revision 1.17  2004/08/03 11:40:35  ivos
      !  Routine tested in 2d and 3d. bugfixes: Added proper type conversions
      !  if KIND and ppm_kind are not the same. Fixed index error in
      !  ppm_psendbuffer.
      !
      !  Revision 1.16  2004/07/26 07:42:45  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.15  2004/06/10 16:20:01  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.14  2004/05/28 11:54:05  walther
      !  Cosmetics.
      !
      !  Revision 1.13  2004/02/18 15:00:32  walther
      !  Renamed the ppm_param_map_update to ppm_param_map_partial.
      !
      !  Revision 1.12  2004/02/05 09:01:21  walther
      !  Added a few more comments.
      !
      !  Revision 1.11  2004/01/26 12:36:52  ivos
      !  Changed ppm_buffer2part from a 2D array (npart,topoid) to a 1D array
      !  (Npart), since it does not depend on the topoid (only one mapping at the
      !  time can be pending since the sendbuffer itself is only 1D as well!).
      !  In map_part_push: removed topoid from the argument list in effect as
      !  it is no longer needed (only 1 mapping can be pending).
      !
      !  Revision 1.10  2004/01/23 17:24:16  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.9  2004/01/23 11:31:22  ivos
      !  Cleanup: (1) updated headers, (2) inserted ppm_error and ppm_write,
      !  (3) added argument checking, (4) added checks after every alloc.
      !
      !  Revision 1.8  2003/12/17 16:52:20  ivos
      !  Changed info on error from 1 to -1.
      !
      !  Revision 1.7  2003/12/16 12:20:45  ivos
      !  Changed alloc_grow to alloc_fit (new mappings can be smaller than 
      !  old ones) and added DEALLOC for ineighsubs at the end of the 
      !  subroutine.
      !
      !  Revision 1.6  2003/12/12 18:02:33  ivos
      !  Changed name of map type to ppm_param_map_update.
      !
      !  Revision 1.5  2003/12/12 16:40:11  ivos
      !  Always uses the internally stored current topology ID and does not 
      !  take topoid as an argument any more.
      !
      !  Revision 1.4  2003/12/11 13:54:31  ivos
      !  Bugfixes: (1) single and double precision arrays were switched in the
      !  #ifdefs when storing the particles to the send buffer. (2) the last 
      !  index in ppm_psendbuffer was corrected from ppm_nneighlist(topoid)+1 
      !  to ppm_ncommseq(topoid)+1.
      !
      !  Revision 1.3  2003/12/09 09:53:30  ivos
      !  Now implements the communication protocol as stored in ppm_icommseq().
      !
      !  Revision 1.2  2003/12/05 11:54:25  ivos
      !  Bugfix. Now compiles.
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_partial_s(xp,Npart,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_partial_d(xp,Npart,info)
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
      INTEGER                 , INTENT(IN   ) :: Npart
      ! list and number of neighboring processors
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                        :: topoid
      INTEGER, DIMENSION(3)          :: ldu
      INTEGER, DIMENSION(:), POINTER :: bcdef
      INTEGER                        :: i,j,k,idom,ipart,nlist1,nlist2
      INTEGER                        :: sendrank,recvrank
      INTEGER                        :: nneighsubs, jdom
      INTEGER                        :: iopt,iset,ibuffer,isonneigh
      INTEGER                        :: recvidx
      CHARACTER(ppm_char)            :: mesg
      REAL(MK)                       :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_partial',t0,info)
      topoid = ppm_topoid
      bcdef => ppm_bcdef(:,topoid)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_partial',  &
     &            'Npart must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

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
      ldu(1) = ppm_nsubs(topoid)
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
      DO i=1,ppm_nsubs(topoid)
          ! if either on local ...
          IF (ppm_subs2proc(i,topoid) .EQ. ppm_rank) THEN
              nneighsubs             = nneighsubs + 1
              ineighsubs(nneighsubs) = i
          ELSE
              ! or on any of my neighbors
              isonneigh = 0
              j = 1
              DO WHILE (isonneigh .EQ. 0 .AND. j .LE. ppm_nneighlist(topoid))
                  IF (ppm_subs2proc(i,topoid).EQ.ppm_ineighlist(j,topoid)) THEN
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
      ldu(1) = ppm_ncommseq(topoid) + 1
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
            sendrank = ppm_subs2proc(idom,topoid)

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
               IF (xp(1,ipart).GE.ppm_min_subs(1,idom,topoid).AND.   &
     &             xp(1,ipart).LE.ppm_max_subs(1,idom,topoid).AND.   &
     &             xp(2,ipart).GE.ppm_min_subs(2,idom,topoid).AND.   & 
     &             xp(2,ipart).LE.ppm_max_subs(2,idom,topoid)) THEN
                   !------------------------------------------------------------
                   !  In the non-periodic case, allow particles that are
                   !  exactly ON an upper EXTERNAL boundary.
                   !------------------------------------------------------------
                   IF((xp(1,ipart).LT.ppm_max_subs(1,idom,topoid) .OR.  &
     &                (ppm_subs_bc(2,idom,topoid).EQ.1           .AND.  &
     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
     &                (xp(2,ipart).LT.ppm_max_subs(2,idom,topoid) .OR.  &
     &                (ppm_subs_bc(4,idom,topoid).EQ.1           .AND.  &
     &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
#elif  __KIND == __DOUBLE_PRECISION
               IF (xp(1,ipart).GE.ppm_min_subd(1,idom,topoid).AND.   &
     &             xp(1,ipart).LE.ppm_max_subd(1,idom,topoid).AND.   &
     &             xp(2,ipart).GE.ppm_min_subd(2,idom,topoid).AND.   & 
     &             xp(2,ipart).LE.ppm_max_subd(2,idom,topoid)) THEN
                   !------------------------------------------------------------
                   !  In the non-periodic case, allow particles that are
                   !  exactly ON an upper EXTERNAL boundary.
                   !------------------------------------------------------------
                   IF((xp(1,ipart).LT.ppm_max_subd(1,idom,topoid) .OR.  &
     &                (ppm_subs_bc(2,idom,topoid).EQ.1           .AND.  &
     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
     &                (xp(2,ipart).LT.ppm_max_subd(2,idom,topoid) .OR.  &
     &                (ppm_subs_bc(4,idom,topoid).EQ.1           .AND.  &
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
            sendrank   = ppm_subs2proc(idom,topoid)
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
               IF (xp(1,ipart).GE.ppm_min_subs(1,idom,topoid).AND.   &
     &             xp(1,ipart).LE.ppm_max_subs(1,idom,topoid).AND.   &
     &             xp(2,ipart).GE.ppm_min_subs(2,idom,topoid).AND.   & 
     &             xp(2,ipart).LE.ppm_max_subs(2,idom,topoid).AND.   &
     &             xp(3,ipart).GE.ppm_min_subs(3,idom,topoid).AND.   & 
     &             xp(3,ipart).LE.ppm_max_subs(3,idom,topoid)) THEN
                   !------------------------------------------------------------
                   !  In the non-periodic case, allow particles that are
                   !  exactly ON an upper EXTERNAL boundary.
                   !------------------------------------------------------------
                   IF((xp(1,ipart).LT.ppm_max_subs(1,idom,topoid) .OR. &
     &                (ppm_subs_bc(2,idom,topoid).EQ.1           .AND. &
     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND. &
     &                (xp(2,ipart).LT.ppm_max_subs(2,idom,topoid) .OR. &
     &                (ppm_subs_bc(4,idom,topoid).EQ.1           .AND. &
     &                bcdef(4).NE. ppm_param_bcdef_periodic))    .AND. &
     &                (xp(3,ipart).LT.ppm_max_subs(3,idom,topoid) .OR. &
     &                (ppm_subs_bc(6,idom,topoid).EQ.1           .AND. &
     &                bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN
#elif  __KIND == __DOUBLE_PRECISION
               IF (xp(1,ipart).GE.ppm_min_subd(1,idom,topoid).AND.   &
     &             xp(1,ipart).LE.ppm_max_subd(1,idom,topoid).AND.   &
     &             xp(2,ipart).GE.ppm_min_subd(2,idom,topoid).AND.   & 
     &             xp(2,ipart).LE.ppm_max_subd(2,idom,topoid).AND.   &
     &             xp(3,ipart).GE.ppm_min_subd(3,idom,topoid).AND.   & 
     &             xp(3,ipart).LE.ppm_max_subd(3,idom,topoid)) THEN
                   !------------------------------------------------------------
                   !  In the non-periodic case, allow particles that are
                   !  exactly ON an upper EXTERNAL boundary.
                   !------------------------------------------------------------
                   IF((xp(1,ipart).LT.ppm_max_subd(1,idom,topoid) .OR. &
     &                (ppm_subs_bc(2,idom,topoid).EQ.1           .AND. &
     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND. &
     &                (xp(2,ipart).LT.ppm_max_subd(2,idom,topoid) .OR. &
     &                (ppm_subs_bc(4,idom,topoid).EQ.1           .AND. &
     &                bcdef(4).NE. ppm_param_bcdef_periodic))    .AND. &
     &                (xp(3,ipart).LT.ppm_max_subd(3,idom,topoid) .OR. &
     &                (ppm_subs_bc(6,idom,topoid).EQ.1           .AND. &
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
      ppm_nsendlist = ppm_ncommseq(topoid)
      ppm_nrecvlist = ppm_ncommseq(topoid)
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
      DO i=1,ppm_ncommseq(topoid)
         !----------------------------------------------------------------------
         !  get next neighbor to send/recv to/from.
         !----------------------------------------------------------------------
         sendrank = ppm_icommseq(i,topoid)
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
      ppm_psendbuffer(ppm_ncommseq(topoid)+1) = npart + 1

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
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_partial_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_partial_d
#endif

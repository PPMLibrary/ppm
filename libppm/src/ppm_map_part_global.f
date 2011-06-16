      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_part_global
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine maps the particles onto the given
      !                 topology using a global mapping (i.e. every
      !                 processor communicates with every other).
      !
      !  Input        : xp(:,:)   (F) the position of the particles
      !                 Npart     (I) the number of particles (on the
      !                               local processor)
      !                 topoid    (I) topology identifier (internal numbering)
      !                               of target
      !
      !  Input/output : 
      !
      !  Output       : info      (I)  return status. 0 on success.
      !
      !  Remarks      : The first part of the buffer contains the on processor
      !                 data.
      !
      !                 The storing of ppm_map_type is not used (yet)
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_part_global.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.24  2006/09/04 18:34:50  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.23  2006/02/06 15:29:10  michaebe
      !  fixed bug 0000016
      !
      !  Revision 1.22  2005/12/01 08:27:19  ivos
      !  For non-periodic cases, particles ON the upper boundary of the
      !  computational domain are now properly accounted for.
      !
      !  Revision 1.21  2005/10/12 11:12:19  ivos
      !  To facilitate non-periodic b.c., particles outside of the domain
      !  are now simply lost.
      !
      !  Revision 1.20  2004/11/11 15:26:17  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.19  2004/10/01 16:09:07  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.18  2004/08/31 12:48:10  ivos
      !  Changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.17  2004/08/30 07:49:18  ivos
      !  Relaxed argument check from Npart.LE.0 to Npart.LT.0.
      !
      !  Revision 1.16  2004/08/03 12:48:28  ivos
      !  Explicitly split loops for ppm_dim.EQ.2 and ppm_dim.EQ.3 cases.
      !  All loops now vectorize on the NEC.
      !
      !  Revision 1.15  2004/08/03 11:39:40  ivos
      !  bugfix: if type of xp and ppm_kind were not the same, the storing
      !  of the particle failed. Fixed by adding proper type conversions.
      !
      !  Revision 1.14  2004/07/26 07:42:45  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.13  2004/06/10 16:20:01  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.12  2004/05/28 11:55:02  walther
      !  Questioned a line in the code.
      !
      !  Revision 1.11  2004/02/20 10:06:49  walther
      !  Bug fix: the ppm_psendbuffer(i+1) should be equal to iset + 1
      !
      !  Revision 1.10  2004/01/26 12:36:51  ivos
      !  Changed ppm_buffer2part from a 2D array (npart,topoid) to a 1D array
      !  (Npart), since it does not depend on the topoid (only one mapping at the
      !  time can be pending since the sendbuffer itself is only 1D as well!).
      !  In map_part_push: removed topoid from the argument list in effect as
      !  it is no longer needed (only 1 mapping can be pending).
      !
      !  Revision 1.9  2004/01/23 17:24:16  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.8  2004/01/23 11:31:22  ivos
      !  Cleanup: (1) updated headers, (2) inserted ppm_error and ppm_write,
      !  (3) added argument checking, (4) added checks after every alloc.
      !
      !  Revision 1.7  2003/12/16 12:23:43  ivos
      !  Replaced alloc_grow with alloc_fit for send buffers (new map can be 
      !  smaller than old ones).
      !
      !  Revision 1.6  2003/12/11 17:30:15  hiebers
      !  removed loop over all processors in serial version
      !
      !  Revision 1.5  2003/12/11 13:55:28  ivos
      !  Bugfix: single and double precision arrays were switched in the #ifdefs
      !  when storing the particles to the send buffer. (removed SIGSEGV).
      !
      !  Revision 1.4  2003/12/10 09:25:58  walther
      !  Updated header and now calling MPI_Abort if a processors fails to
      !  map the particles.
      !
      !  Revision 1.3  2003/12/05 09:45:08  ivos
      !  Condensed building of send/recv buffers into a single loop using 
      !  cpp #if for the different precision types.
      !
      !  Revision 1.2  2003/12/05 09:18:14  ivos
      !  now using the internally stored ppm_dim instead of taking it as an
      !  argument.
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
      SUBROUTINE ppm_map_part_global_s(xp,Npart,topoid,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_global_d(xp,Npart,topoid,info)
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
      USE ppm_module_check_topoid
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
      INTEGER                 , INTENT(IN   ) :: Npart,topoid
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)          :: ldu
      INTEGER, DIMENSION(:), POINTER :: bcdef
      INTEGER                        :: i,j,k,idom,ipart,nlist1,nlist2
      INTEGER                        :: sendrank,recvrank
      INTEGER                        :: iopt,iset,ibuffer
      CHARACTER(ppm_char)            :: mesg
      REAL(MK)                       :: t0
      LOGICAL                        :: valid
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_global',t0,info)
      bcdef => ppm_bcdef(:,topoid)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_global',  &
     &            'Npart must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_check_topoid(ppm_param_id_internal,topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_global',  &
     &            'topoid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls (not used yet)
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_global

      !-------------------------------------------------------------------------
      !  Alloc memory for particle lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = Npart
      CALL ppm_alloc(ilist1,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_global',     &
     &        'particle list 1 ILIST1',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_global',     &
     &        'particle list 2 ILIST2',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the pointer to the buffer; for the global map we
      !  need entries for each processor, thus ldu(1) = ppm_nproc
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nproc + 1
      CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_global',     &
     &        'particle send buffer PPM_PSENDBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      iopt   = ppm_param_alloc_fit
      ldu(1) = Npart
      CALL ppm_alloc(part2proc,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_global',     &
     &        'particles-to-processor map PART2PROC',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = Npart
      CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_global',     &
     &        'buffer-to-particles map PPM_BUFFER2PART',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the particle list
      !-------------------------------------------------------------------------
      nlist1 = 0
      DO ipart=1,Npart
         nlist1         = nlist1 + 1
         ilist1(nlist1) = ipart
      ENDDO

      !-------------------------------------------------------------------------
      !  Then of these exclude the particles that on the current processor
      !-------------------------------------------------------------------------
      IF (ppm_dim.EQ.2) THEN
         !----------------------------------------------------------------------
         !  Loop over the subdomains 
         !----------------------------------------------------------------------
!        DO idom=1,ppm_nsubs(topoid)
         DO idom=ppm_nsubs(topoid),1,-1
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
         !  nlist2 as fast as possible)
         !----------------------------------------------------------------------
!        DO idom=1,ppm_nsubs(topoid)
         DO idom=ppm_nsubs(topoid),1,-1
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
      !  Here we could check that we sold all the particles, but if we use 
      !  Dirichlet BCs we might just want to loose the particles outside of the 
      !  domain. So just skip them in the map. This is also consistent with the 
      !  map_partial where particles are also lost (there are no periodic sub 
      !  images in this case).
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Store the number of buffer entries (this is the first)
      !-------------------------------------------------------------------------
      ppm_buffer_set = 1

      !-------------------------------------------------------------------------
      !  Allocate memory for the field registers that holds the dimension and 
      !  type of the data
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_buffer_set
      CALL ppm_alloc(ppm_buffer_dim ,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_global',     &
     &        'buffer dimensions PPM_BUFFER_DIM',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_global',     &
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
      !  (Re)allocate memory for the buffer (need not save the contents of the
      !  buffer since this is a global map).
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
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_global',     &
     &        'global send buffer PPM_SENDBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the sendlist
      !-------------------------------------------------------------------------
      ppm_nsendlist = ppm_nproc 
      ppm_nrecvlist = ppm_nproc 
      ldu(1)        = ppm_nsendlist
      CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_global',     &
     &        'send list PPM_ISENDLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_global',     &
     &        'receive list PPM_IRECVLIST',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the particle lists 
      !-------------------------------------------------------------------------
      nlist1 = npart
      DO ipart=1,Npart
         ilist1(ipart) = ipart
      ENDDO

      !-------------------------------------------------------------------------
      !  loop over all processors, starting with the processor itself 
      !-------------------------------------------------------------------------
      sendrank           = ppm_rank - 1
      recvrank           = ppm_rank + 1
      ppm_psendbuffer(1) = 1
      ppm_nsendlist      = 0
      ppm_nrecvlist      = 0
      iset               = 0
      ibuffer            = 0

!SEH: loop over all processors only in parallel version needed
#if defined __MPI
      DO i=1,ppm_nproc
         !----------------------------------------------------------------------
         !  compute the next processor
         !----------------------------------------------------------------------
         sendrank = sendrank + 1
         IF (sendrank.GT.ppm_nproc-1) sendrank = sendrank - ppm_nproc 
         recvrank = recvrank - 1
         IF (recvrank.LT.          0) recvrank = recvrank + ppm_nproc 

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
         !  Initialize the buffer count
         !----------------------------------------------------------------------
         nlist2 = 0

         IF (ppm_dim .EQ. 2) THEN
            DO j=1,nlist1
               ipart = ilist1(j)
               IF (part2proc(ipart).EQ.sendrank) THEN  
                  !-------------------------------------------------------------
                  !  increment the buffer counter
                  !-------------------------------------------------------------
                  iset = iset + 1
  
                  !-------------------------------------------------------------
                  !  Store the id of the particle
                  !-------------------------------------------------------------
                  ppm_buffer2part(iset) = ipart
    
                  !-------------------------------------------------------------
                  !  Store the particle
                  !-------------------------------------------------------------
                  IF (ppm_kind .EQ. ppm_kind_double) THEN
#if    __KIND == __SINGLE_PRECISION
                     ibuffer                  = ibuffer + 1
                     ppm_sendbufferd(ibuffer) = REAL(xp(1,ipart),   &
     &                   ppm_kind_double)
                     ibuffer                  = ibuffer + 1
                     ppm_sendbufferd(ibuffer) = REAL(xp(2,ipart),   &
     &                   ppm_kind_double)
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
     &                   ppm_kind_single)
                     ibuffer                  = ibuffer + 1
                     ppm_sendbuffers(ibuffer) = REAL(xp(2,ipart),   &
     &                   ppm_kind_single)
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
                  !-------------------------------------------------------------
                  !  increment the buffer counter
                  !-------------------------------------------------------------
                  iset = iset + 1
  
                  !-------------------------------------------------------------
                  !  Store the id of the particle
                  !-------------------------------------------------------------
                  ppm_buffer2part(iset) = ipart
    
                  !-------------------------------------------------------------
                  !  Store the particle
                  !-------------------------------------------------------------
                  IF (ppm_kind .EQ. ppm_kind_double) THEN
#if    __KIND == __SINGLE_PRECISION
                     ibuffer                  = ibuffer + 1
                     ppm_sendbufferd(ibuffer) = REAL(xp(1,ipart),   &
     &                   ppm_kind_double)
                     ibuffer                  = ibuffer + 1
                     ppm_sendbufferd(ibuffer) = REAL(xp(2,ipart),   &
     &                   ppm_kind_double)
                     ibuffer                  = ibuffer + 1
                     ppm_sendbufferd(ibuffer) = REAL(xp(3,ipart),   &
     &                   ppm_kind_double)
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
     &                   ppm_kind_single)
                     ibuffer                  = ibuffer + 1
                     ppm_sendbuffers(ibuffer) = REAL(xp(2,ipart),   &
     &                   ppm_kind_single)
                     ibuffer                  = ibuffer + 1
                     ppm_sendbuffers(ibuffer) = REAL(xp(3,ipart),   &
     &                   ppm_kind_single)
#endif
                  ENDIF
               ELSE
                  nlist2         = nlist2 + 1
                  ilist2(nlist2) = ipart
               ENDIF
            ENDDO
         ENDIF

         !----------------------------------------------------------------------
         !  Update the buffer pointer
         !----------------------------------------------------------------------
         ppm_psendbuffer(i+1) = iset + 1 

         !----------------------------------------------------------------------
         !  Swap the lists
         !----------------------------------------------------------------------
         nlist1 = nlist2
         DO j=1,nlist1
            ilist1(j) = ilist2(j)
         ENDDO
      ENDDO
#endif
      !-------------------------------------------------------------------------
      !  Everything has to go (should not be necessary!)
      !  the line above: 
      !    
      !      ppm_psendbuffer(i+1) = iset + 1 
      !    
      !  should take care of this
      !-------------------------------------------------------------------------
      ppm_psendbuffer(ppm_nproc+1) = npart + 1

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
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_global',     &
     &        'particle list 1 ILIST1',__LINE__,info)
      ENDIF
      CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_global',     &
     &        'particle list 2 ILIST2',__LINE__,info)
      ENDIF
      CALL ppm_alloc(part2proc,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_global',     &
     &        'particle-to-processor map PART2PROC',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_global',t0,info)
      RETURN
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_global_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_global_d
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_map_part_ghost_get
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine maps/adds the ghost particles on the 
      !                 current topology. This routine is similar to the 
      !                 partial mapping routine (ppm_map_part_partial) in the 
      !                 sense that the ghost particles are assumed to be 
      !                 located on neighbouring processors only, and thus only
      !                 require a nearest neighbour communication. 
      ! 
      !  Input        : xp(:,:)      (F) : the position of the particles
      !                 lda          (I) : leading dimension of xp
      !                 Npart        (I) : the number of particles (on the
      !                                    local processor)
      !                 isymm        (I) : indicator for the use of symmetry 
      !                                    isymm > 0 use symmetry
      !                                    isymm = 0 do not use symmetry
      !                 ghostsize    (F) : the size of the ghost layer
      !                                    
      !  Input/output : info         (I) : return status, 0 on success
      !
      !  Remarks      : This routine SHOULD be efficient - since it will be 
      !                 called frequently. One way of improving (?) the
      !                 performance is to use the cell lists to find the 
      !                 potential ghosts.
      !
      !                 This routine should be split in several routines !
      !
      !                 The ghosts are found in three steps: 
      !
      !                    1) consider ghosts within the local processor
      !                    2) consider ghosts on neighbouring processors
      !                    1) consider ghosts on periodic images of neighbouring
      !                       processors
      !                  
      !                 Comments:
      !                  
      !                 1) local ghosts comes in two types: ghost that exists
      !                    because two sub domains touch, and ghosts across a 
      !                    periodic boundary. The first is automatically handled
      !                    by the particles - and no copy or/and send is 
      !                    required. The second ghosts could be handled during 
      !                    the calculation of the interactions, but this would 
      !                    require a check for periodicity and no explicit 
      !                    ghosts from the processor itself. However, it seems
      !                    more natural to have ghosts - irrespectively of their
      !                    origin, so ghosts from periodicity (on the same 
      !                    processor) are copied here.
      !
      !                 2) is standard procedure
      ! 
      !                 3) is handled by copying particles that are adjacent to
      !                    faces on a physical boundary to their image position.
      !
      !                 Now, item 2) and 3) can be treated within the same logic
      !                 whereas 1) require a bit of thought: to keep the source
      !                 compact, we could loop over all processors including
      !                 the local one and check for ghosts - the problem with 
      !                 this procedure is, that we check for ghosts by comparing
      !                 the location of particles within the extended subs 
      !                 boundaries (and not their true ghost layer - which 
      !                 would result in more IFs than necessary). However, 
      !                 because of this, looping over the processor itself 
      !                 would find ghosts that are not really ghost - the only 
      !                 ghosts that should be found in this step are those due 
      !                 to periodicity.  The solution is to store the ghosts as
      !                 ighost() with nghost denoting the total number of 
      !                 ghosts (excluding those due to periodicity) and 
      !                 nghostplus the total number of ghosts. During the loop 
      !                 over the processor itself we therefore only need to 
      !                 loop from nghost+1,nghostplus.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_part_ghost_get.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.17  2006/10/10 20:51:47  walther
      !  Now added the storing of the ppm_ghost_offsets/d.
      !  Even stored if not needed - but would require too
      !  many IF statements to avoid !?
      !
      !  Revision 1.16  2006/09/04 18:34:50  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.15  2006/02/03 09:41:26  ivos
      !  Added the PRELIMINARY ghost_put functionality. Still needs clean-up,
      !  but should work.
      !
      !  Revision 1.14  2005/08/24 15:32:45  hiebers
      !  BUGFIX: nlist2 is also initialized outside loop for the case that
      !  ppm_nsublist = 0
      !
      !  Revision 1.13  2004/11/11 15:26:17  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.12  2004/10/01 16:09:06  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.11  2004/08/27 11:40:20  walther
      !  Bug fix in ghosts at right,upper,top boundaries.
      !
      !  Revision 1.10  2004/08/24 12:09:00  walther
      !  Bug fix in allocation of ighost() for periodic systems.
      !
      !  Revision 1.9  2004/08/03 12:46:38  ivos
      !  bugfix: if KIND and ppm_kind were not the same, things went wrong.
      !  Fixed by adding proper IFs and type conversions.
      !
      !  Revision 1.8  2004/07/26 07:42:44  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.7  2004/07/15 14:16:57  walther
      !  Bug fix: a ghost could incorrectly be send more than once to a 
      !  processor; introducing the lghost(:) array fixed this problem.
      !
      !  Revision 1.6  2004/07/07 13:09:35  walther
      !  Bug fix: particles in a sub require GE and LT.
      !
      !  Revision 1.5  2004/07/01 13:15:23  walther
      !  Bug fix: allocation of ppm_sendbuffer and ppm_buffer2part.
      !
      !  Revision 1.4  2004/06/03 16:07:22  ivos
      !  Removed debug PRINT statement.
      !
      !  Revision 1.3  2004/06/03 08:49:39  walther
      !  Bug fix: now reallocating the ppm_buffer2part and ppm_sendbuffers/d.
      !  Also skipping the do-loops in step 1 if if we do not have periodicity.
      !
      !  Revision 1.2  2004/05/28 10:32:56  walther
      !  First functional release - without symmetry.
      !
      !  Revision 1.1  2004/02/19 15:56:35  walther
      !  Initial implementation (not complete!).
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_get_s(xp,lda,Npart,isymm,ghostsize,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_get_d(xp,lda,Npart,isymm,ghostsize,info)
#endif 

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION  | __KIND_AUX == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      REAL(MK)                , INTENT(IN   ) :: ghostsize
      INTEGER                 , INTENT(IN   ) :: Npart,lda
      INTEGER                 , INTENT(IN   ) :: isymm
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: i,j,k,topoid,isub
      INTEGER               :: nlist1,nlist2,nghost,nghostplus
      INTEGER               :: ipart,sendrank,recvrank
      INTEGER               :: iopt,iset,ibuffer
      REAL(MK), DIMENSION(:,:), POINTER :: xt  ! position of potential ghosts
      REAL(MK), DIMENSION(:,:), POINTER :: xt_offset ! offset of pot. ghosts
      REAL(MK)              :: xminf,yminf,zminf ! full domain
      REAL(MK)              :: xmaxf,ymaxf,zmaxf ! full domain
      REAL(MK)              :: xmini,ymini,zmini ! inner domain
      REAL(MK)              :: xmaxi,ymaxi,zmaxi ! inner domain
      REAL(MK), DIMENSION(ppm_dim) :: len_phys
      REAL(MK)              :: t0
      ! number of periodic directions: between 0 and ppm_dim
      INTEGER               :: iperiodic
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_ghost_get',t0,info)
      topoid = ppm_topoid

      !-------------------------------------------------------------------------
      !  Compute the size of the computational box on this topology
      !-------------------------------------------------------------------------
#if __KIND == __DOUBLE_PRECISION
      len_phys(:) = ppm_max_physd(:,topoid) - ppm_min_physd(:,topoid)
#else
      len_phys(:) = ppm_max_physs(:,topoid) - ppm_min_physs(:,topoid)
#endif

      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls 
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_ghost_get

      !-------------------------------------------------------------------------
      !  Allocate memory for the list of particle on the local processor that
      !  may be ghosts on other processors
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = Npart
      CALL ppm_alloc(ilist1,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &        'list1',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &        'list2',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ighost,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &        'ighost',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for their position
      !-------------------------------------------------------------------------
      ldu(1) = ppm_dim
      ldu(2) = Npart
      CALL ppm_alloc(xt,ldu,iopt,info)
      CALL ppm_alloc(xt_offset,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &        'xt',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  List ilist1() holds the particles we are currently considering 
      !  List ilist2() holds the particles that have not yet been associated 
      !  with a sub, and ghost() holds the ghosts
      !-------------------------------------------------------------------------
      DO i=1,Npart
         ilist1(i) = i
      ENDDO
      nlist1  = Npart
      nlist2  = 0
      nghost  = 0

      !-------------------------------------------------------------------------
      !  Fill the list with particles that are within the reach of the ghost
      !  regions of other subs. We do that by looping over the subs belonging
      !  to this processor and checking if the particles are well within the
      !  sub. If this is the case, the particle will never be a ghost 
      !-------------------------------------------------------------------------
      DO k=1,ppm_nsublist(topoid)
         !----------------------------------------------------------------------
         !  Initialize the second list counter to zero
         !----------------------------------------------------------------------
         nlist2 = 0

         !----------------------------------------------------------------------
         !  Get the (global) id of the sub
         !----------------------------------------------------------------------
         isub = ppm_isublist(k,topoid)

         !----------------------------------------------------------------------
         !  Store the full extend of the sub in ?minf and ?maxf
         !----------------------------------------------------------------------
#if __KIND == __DOUBLE_PRECISION
         xminf = ppm_min_subd(1,isub,topoid)
         xmaxf = ppm_max_subd(1,isub,topoid)

         yminf = ppm_min_subd(2,isub,topoid)
         ymaxf = ppm_max_subd(2,isub,topoid)

         IF (ppm_dim.EQ.3) THEN
            zminf = ppm_min_subd(3,isub,topoid)
            zmaxf = ppm_max_subd(3,isub,topoid)
         ENDIF 
#else
         xminf = ppm_min_subs(1,isub,topoid)
         xmaxf = ppm_max_subs(1,isub,topoid)

         yminf = ppm_min_subs(2,isub,topoid)
         ymaxf = ppm_max_subs(2,isub,topoid)
         IF (ppm_dim.EQ.3) THEN
            zminf = ppm_min_subs(3,isub,topoid)
            zmaxf = ppm_max_subs(3,isub,topoid)
         ENDIF
#endif

         !----------------------------------------------------------------------
         !  Compute the size of the inner region
         !----------------------------------------------------------------------
         IF (isymm.GT.0) THEN
            !-------------------------------------------------------------------
            !  if we use symmetry the upper/right part of our sub will have 
            !  ghosts and therefore the particle at the upper/right part of the 
            !  sub cannot be ghosts on other processors. Thus the ghosts must 
            !  be found at the lower/left of the sub
            !-------------------------------------------------------------------
            xmini = xminf + ghostsize
            xmaxi = xmaxf

            ymini = yminf + ghostsize
            ymaxi = ymaxf
 
            IF (ppm_dim.EQ.3) THEN
               zmini = zminf + ghostsize
               zmaxi = zmaxf
            ENDIF 
         ELSE
            !-------------------------------------------------------------------
            !  if we do not use symmetry, the particles along the entire 
            !  boundary of the sub will be ghosts on other processors
            !-------------------------------------------------------------------
            xmini = xminf + ghostsize
            xmaxi = xmaxf - ghostsize

            ymini = yminf + ghostsize
            ymaxi = ymaxf - ghostsize

            IF (ppm_dim.EQ.3) THEN
               zmini = zminf + ghostsize
               zmaxi = zmaxf - ghostsize
            ENDIF 
         ENDIF 

         !----------------------------------------------------------------------
         !  loop over the remaining particles
         !----------------------------------------------------------------------
         IF (ppm_dim.EQ.2) THEN
            DO j=1,nlist1
               !----------------------------------------------------------------
               !  get the particle index
               !----------------------------------------------------------------
               ipart = ilist1(j)
   
               !----------------------------------------------------------------
               !  check if the particles belongs to this sub
               !----------------------------------------------------------------
               IF (xp(1,ipart).GE.xminf.AND.xp(1,ipart).LT.xmaxf.AND. &
     &             xp(2,ipart).GE.yminf.AND.xp(2,ipart).LT.ymaxf) THEN
                  !-------------------------------------------------------------
                  !  if yes, check if the particles are in the ghost layer 
                  !-------------------------------------------------------------
                  IF (xp(1,ipart).LT.xmini.OR.xp(1,ipart).GT.xmaxi.OR. &
     &                xp(2,ipart).LT.ymini.OR.xp(2,ipart).GT.ymaxi) THEN
                     !----------------------------------------------------------
                     !  if yes the particle will be a ghost somewhere
                     !----------------------------------------------------------
                     nghost         = nghost + 1   
                     ighost(nghost) = ipart
                     xt(1,nghost)   = xp(1,ipart)
                     xt(2,nghost)   = xp(2,ipart)
                     xt_offset(1,nghost) = 0.0_MK
                     xt_offset(2,nghost) = 0.0_MK
                  ENDIF
               ELSE    
                  !-------------------------------------------------------------
                  !  If not on this sub we need to consider it further
                  !-------------------------------------------------------------
                  nlist2         = nlist2 + 1
                  ilist2(nlist2) = ipart
               ENDIF
            ENDDO
         ELSE
            DO j=1,nlist1
               !----------------------------------------------------------------
               !  get the particle index
               !----------------------------------------------------------------
               ipart = ilist1(j)
   
               !----------------------------------------------------------------
               !  check if the particles belongs to this sub
               !----------------------------------------------------------------
               IF (xp(1,ipart).GE.xminf.AND.xp(1,ipart).LT.xmaxf.AND. &
     &             xp(2,ipart).GE.yminf.AND.xp(2,ipart).LT.ymaxf.AND. &
     &             xp(3,ipart).GE.zminf.AND.xp(3,ipart).LT.zmaxf) THEN
                  !-------------------------------------------------------------
                  !  if yes, check if the particles is in the ghost layer 
                  !-------------------------------------------------------------
                  IF (xp(1,ipart).LT.xmini.OR.xp(1,ipart).GT.xmaxi.OR. &
     &                xp(2,ipart).LT.ymini.OR.xp(2,ipart).GT.ymaxi.OR. &
     &                xp(3,ipart).LT.zmini.OR.xp(3,ipart).GT.zmaxi) THEN
                     !----------------------------------------------------------
                     !  if yes the particle will be a ghost somewhere
                     !----------------------------------------------------------
                     nghost         = nghost + 1   
                     ighost(nghost) = ipart
                     xt(1,nghost)   = xp(1,ipart)
                     xt(2,nghost)   = xp(2,ipart)
                     xt(3,nghost)   = xp(3,ipart)
                     xt_offset(1,nghost) = 0.0_MK
                     xt_offset(2,nghost) = 0.0_MK
                     xt_offset(3,nghost) = 0.0_MK
                  ENDIF
               ELSE    
                  !-------------------------------------------------------------
                  !  If not on this sub we need to consider it further
                  !-------------------------------------------------------------
                  nlist2         = nlist2 + 1
                  ilist2(nlist2) = ipart
               ENDIF
            ENDDO
         ENDIF 

         !----------------------------------------------------------------------
         !  swap the lists
         !----------------------------------------------------------------------
         DO j=1,nlist2
            ilist1(j) = ilist2(j)
         ENDDO
         nlist1 = nlist2

      ENDDO ! end of subs on local processor



      !-------------------------------------------------------------------------
      !  At the end the nlist2 should be zero
      !-------------------------------------------------------------------------
      IF (nlist2.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_part_unass,'ppm_map_part_ghost_get',     &
     &       'nlist2 > 0',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the total number of ghosts incl. those due to periodicity
      !-------------------------------------------------------------------------
      nghostplus = nghost

      !-------------------------------------------------------------------------
      !  Ok, we now have a list of potential ghosts. From these we extract/add
      !  their periodic images (if any). So, first we check for periodicity
      !-------------------------------------------------------------------------
      iperiodic = 0
      DO k=1,ppm_dim
         IF (ppm_bcdef(2*k-1,topoid).EQ.ppm_param_bcdef_periodic) THEN
            iperiodic = iperiodic + 1
         ENDIF 
      ENDDO

      !-------------------------------------------------------------------------
      !  If we have periodicity, create the periodic ghosts
      !-------------------------------------------------------------------------
      IF (iperiodic.GT.0) THEN
         !----------------------------------------------------------------------
         !  handle periodicity in x
         !----------------------------------------------------------------------
         IF (ppm_bcdef(1,topoid).EQ.ppm_param_bcdef_periodic) THEN
            !-------------------------------------------------------------------
            !  (Re)allocate memory for the periodic ghosts
            !-------------------------------------------------------------------
            iopt   = ppm_param_alloc_grow_preserve
            ldu(1) = ppm_dim
            ldu(2) = 2*nghostplus
            CALL ppm_alloc(xt,ldu,iopt,info) 
            CALL ppm_alloc(xt_offset,ldu,iopt,info) 
            IF (info.NE.0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &             'xt',__LINE__,info)
               GOTO 9999
            ENDIF

            ldu(1) = ldu(2)
            CALL ppm_alloc(ighost,ldu,iopt,info) 
            IF (info.NE.0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &              'ighost',__LINE__,info)
                GOTO 9999
            ENDIF

            !-------------------------------------------------------------------
            !  copy periodic ghosts in the x-direction
            !-------------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
            xminf = ppm_min_physs(1,topoid) 
            xmini = ppm_min_physs(1,topoid) + ghostsize
#else
            xminf = ppm_min_physd(1,topoid) 
            xmini = ppm_min_physd(1,topoid) + ghostsize
#endif
            k     = nghostplus
            DO i=1,nghostplus
               !----------------------------------------------------------------
               !  first those at the west boundary 
               !----------------------------------------------------------------
               IF (xt(1,i).GE.xminf.AND.xt(1,i).LT.xmini) THEN
                  k         = k + 1
                  ighost(k) = ighost(i)
                  xt(1,k)   = xt(1,i) + len_phys(1)
                  xt(2,k)   = xt(2,i)
                  xt_offset(1,k) = len_phys(1)
                  xt_offset(2,k) = 0.0_MK
                  IF (ppm_dim.EQ.3) THEN
                     xt(3,k)   = xt(3,i)
                     xt_offset(3,k) = 0.0_MK
                  ENDIF 
               ENDIF
            ENDDO
            IF (isymm.EQ.0) THEN
               !----------------------------------------------------------------
               !  then the east bc, but only if we are not using symmetry
               !----------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
               xmaxf = ppm_max_physs(1,topoid) 
               xmaxi = ppm_max_physs(1,topoid) - ghostsize
#else
               xmaxf = ppm_max_physd(1,topoid) 
               xmaxi = ppm_max_physd(1,topoid) - ghostsize
#endif
               DO i=1,nghostplus
                  IF  (xt(1,i).GT.xmaxi.AND.xt(1,i).LT.xmaxf) THEN
                     k         = k + 1
                     ighost(k) = ighost(i)
                     xt(1,k)   = xt(1,i) - len_phys(1)
                     xt(2,k)   = xt(2,i)
                     xt_offset(1,k) = -len_phys(1)
                     xt_offset(2,k) = 0.0_MK
                     IF (ppm_dim.EQ.3) THEN
                        xt(3,k)   = xt(3,i)
                        xt_offset(3,k) = 0.0_MK
                     ENDIF 
                  ENDIF
               ENDDO
            ENDIF 

            !-------------------------------------------------------------------
            !  update the ghost counter
            !-------------------------------------------------------------------
            nghostplus = k
         ENDIF ! of periodicity in x 

         !----------------------------------------------------------------------
         !  handle periodicity in y
         !----------------------------------------------------------------------
         IF (ppm_bcdef(3,topoid).EQ.ppm_param_bcdef_periodic) THEN
            !-------------------------------------------------------------------
            !  (Re)allocate memory for the periodic ghosts
            !-------------------------------------------------------------------
            iopt   = ppm_param_alloc_grow_preserve
            ldu(1) = ppm_dim
            ldu(2) = 2*nghostplus
            CALL ppm_alloc(xt,ldu,iopt,info) 
            CALL ppm_alloc(xt_offset,ldu,iopt,info) 
            IF (info.NE.0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &             'xt',__LINE__,info)
               GOTO 9999
            ENDIF

            ldu(1) = ldu(2)
            CALL ppm_alloc(ighost,ldu,iopt,info) 
            IF (info.NE.0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &              'ighost',__LINE__,info)
                GOTO 9999
            ENDIF

            !-------------------------------------------------------------------
            !  copy periodic ghosts in the y-direction
            !-------------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
            yminf = ppm_min_physs(2,topoid) 
            ymini = ppm_min_physs(2,topoid) + ghostsize
#else
            yminf = ppm_min_physd(2,topoid) 
            ymini = ppm_min_physd(2,topoid) + ghostsize
#endif
            k     = nghostplus
            DO i=1,nghostplus
               !----------------------------------------------------------------
               !  first those at the south boundary 
               !----------------------------------------------------------------
               IF (xt(2,i).GE.yminf.AND.xt(2,i).LT.ymini) THEN
                  k         = k + 1
                  ighost(k) = ighost(i)
                  xt(1,k)   = xt(1,i) 
                  xt(2,k)   = xt(2,i) + len_phys(2)
                  xt_offset(1,k) = 0.0_MK
                  xt_offset(2,k) = len_phys(2)
                  IF (ppm_dim.EQ.3) THEN
                     xt(3,k)   = xt(3,i)
                     xt_offset(3,k) = 0.0_MK
                  ENDIF 
               ENDIF
            ENDDO
            IF (isymm.EQ.0) THEN
               !----------------------------------------------------------------
               !  then the north bc, but only if we are not using symmetry
               !----------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
               ymaxf = ppm_max_physs(2,topoid) 
               ymaxi = ppm_max_physs(2,topoid) - ghostsize
#else
               ymaxf = ppm_max_physd(2,topoid) 
               ymaxi = ppm_max_physd(2,topoid) - ghostsize
#endif
               DO i=1,nghostplus
                  IF  (xt(2,i).GT.ymaxi.AND.xt(2,i).LT.ymaxf) THEN
                     k         = k + 1
                     ighost(k) = ighost(i)
                     xt(1,k)   = xt(1,i)
                     xt(2,k)   = xt(2,i) - len_phys(2)
                     xt_offset(1,k) = 0.0_MK
                     xt_offset(2,k) = -len_phys(2)
                     IF (ppm_dim.EQ.3) THEN
                        xt(3,k)   = xt(3,i)
                        xt_offset(3,k) = 0.0_MK
                     ENDIF 
                  ENDIF
               ENDDO
            ENDIF 

            !-------------------------------------------------------------------
            !  update the ghost counter
            !-------------------------------------------------------------------
            nghostplus = k
         ENDIF ! of periodicity in y 

         !----------------------------------------------------------------------
         !  handle periodicity in z (if 3D)
         !----------------------------------------------------------------------
         IF (ppm_dim.EQ.3) THEN
            !-------------------------------------------------------------------
            !  yes, we split the if in two, since we do not know in what order
            !  the compiler will check and ppm_bcdef will only be allocated to
            !  four (4) in 2D
            !-------------------------------------------------------------------
            IF (ppm_bcdef(5,topoid).EQ.ppm_param_bcdef_periodic) THEN
               !----------------------------------------------------------------
               !  (Re)allocate memory for the periodic ghosts
               !----------------------------------------------------------------
               iopt   = ppm_param_alloc_grow_preserve
               ldu(1) = ppm_dim
               ldu(2) = 2*nghostplus
               CALL ppm_alloc(xt,ldu,iopt,info) 
               CALL ppm_alloc(xt_offset,ldu,iopt,info) 
               IF (info.NE.0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &                'xt',__LINE__,info)
                  GOTO 9999
               ENDIF

               ldu(1) = ldu(2)
               CALL ppm_alloc(ighost,ldu,iopt,info) 
               IF (info.NE.0) THEN
                   info = ppm_error_fatal
                   CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &                 'ighost',__LINE__,info)
                   GOTO 9999
               ENDIF

               !----------------------------------------------------------------
               !  copy periodic ghosts in the z-direction
               !----------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
               zminf = ppm_min_physs(3,topoid) 
               zmini = ppm_min_physs(3,topoid) + ghostsize
#else
               zminf = ppm_min_physd(3,topoid) 
               zmini = ppm_min_physd(3,topoid) + ghostsize
#endif
               k     = nghostplus 
               DO i=1,nghostplus
                  !-------------------------------------------------------------
                  !  first those at the south boundary 
                  !-------------------------------------------------------------
                  IF (xt(3,i).GE.zminf.AND.xt(3,i).LT.zmini) THEN
                     k         = k + 1
                     ighost(k) = ighost(i)
                     xt(1,k)   = xt(1,i) 
                     xt(2,k)   = xt(2,i)
                     xt(3,k)   = xt(3,i) + len_phys(3)

                     xt_offset(1,k) = 0.0_MK
                     xt_offset(2,k) = 0.0_MK
                     xt_offset(3,k) = len_phys(3)
                  ENDIF
               ENDDO
               IF (isymm.EQ.0) THEN
                  !-------------------------------------------------------------
                  !  then the north bc, but only if we are not using symmetry
                  !-------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
                  zmaxf = ppm_max_physs(3,topoid) 
                  zmaxi = ppm_max_physs(3,topoid) - ghostsize
#else
                  zmaxf = ppm_max_physd(3,topoid) 
                  zmaxi = ppm_max_physd(3,topoid) - ghostsize
#endif
                  DO i=1,nghostplus
                     IF  (xt(3,i).GT.zmaxi.AND.xt(3,i).LT.zmaxf) THEN
                        k         = k + 1
                        ighost(k) = ighost(i)
                        xt(1,k)   = xt(1,i)
                        xt(2,k)   = xt(2,i) 
                        xt(3,k)   = xt(3,i) - len_phys(3)
                     ENDIF
                  ENDDO
               ENDIF 

               !----------------------------------------------------------------
               !  update the ghost counter
               !----------------------------------------------------------------
               nghostplus = k
            ENDIF ! of periodicity in y 
         ENDIF ! of 3D
      ENDIF ! of periodicity at all/any direction

      !-------------------------------------------------------------------------
      !  Ok, now we have a shorter list of potential ghosts to search 
      !  allocate memory for the lghosts 
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = nghostplus
      CALL ppm_alloc(lghost,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &        'logical list: lghost',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the global send/recv lists
      !-------------------------------------------------------------------------
      ldu(1) = ppm_ncommseq(topoid)
      CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &        'global send rank list PPM_ISENDLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &        'global recv rank list PPM_IRECVLIST',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the pointers to the particles that will be send
      !  and received from and to a processor. The ppm_recvbuffer is NOT used
      !  in this routine, but initialized from the precv() array in the routine
      !  ppm_map_part_send().
      !-------------------------------------------------------------------------
      ldu(1) = ppm_ncommseq(topoid) + 1
      CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &        'global send buffer pointer PPM_PSENDBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Step 1:
      !  First we find the ghosts due to periodicity on the local processor
      !-------------------------------------------------------------------------
      ppm_psendbuffer(1) = 1
      ppm_nsendlist      = 0
      ppm_nrecvlist      = 0
      iset               = 0
      ibuffer            = 0
      k                  = 1

      !-------------------------------------------------------------------------
      !  Get the rank of the processor
      !-------------------------------------------------------------------------
      sendrank = ppm_icommseq(k,topoid) ! should be ppm_rank !
      recvrank = sendrank

      !-------------------------------------------------------------------------
      !  Store the processor to which we will send to
      !-------------------------------------------------------------------------
      ppm_nsendlist                = ppm_nsendlist + 1
      ppm_isendlist(ppm_nsendlist) = sendrank

      !-------------------------------------------------------------------------
      !  Store the processor to which we will recv from
      !-------------------------------------------------------------------------
      ppm_nrecvlist                = ppm_nrecvlist + 1
      ppm_irecvlist(ppm_nrecvlist) = recvrank

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
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &        'buffer dimensions PPM_BUFFER_DIM',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &        'buffer types PPM_BUFFER_TYPE',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the dimension and type
      !-------------------------------------------------------------------------
      ppm_buffer_dim(ppm_buffer_set)  = ppm_dim
#if    __KIND == __SINGLE_PRECISION
      ppm_buffer_type(ppm_buffer_set) = ppm_kind_single
#else
      ppm_buffer_type(ppm_buffer_set) = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Well we only do the stuff for real if we have periodicity !
      !  (Re)allocate memory for the send buffer. The required size of these 
      !  arrays is not easy to compute. The first step (step 1) require at
      !  most ppm_nsublist(topoid)*(nghostplus - nghost)*ppm_dim. The arrays
      !  are resized further below during step 2 and 3.
      !-------------------------------------------------------------------------
      IF (iperiodic.GT.0) THEN
         !----------------------------------------------------------------------
         !  Well we can grow or fit the arrays - a matter of taste 
         !----------------------------------------------------------------------
         iopt   = ppm_param_alloc_grow
         ldu(1) = ppm_dim*(nghostplus - nghost)*ppm_nsublist(topoid)
         !----------------------------------------------------------------------
         !  First allocate the sendbuffer
         !----------------------------------------------------------------------
         IF (ppm_kind.EQ.ppm_kind_double) THEN
            CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
            CALL ppm_alloc(ppm_ghost_offsetd,ldu,iopt,info)
         ELSE
            CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
            CALL ppm_alloc(ppm_ghost_offsets,ldu,iopt,info)
         ENDIF
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &           'global send buffer PPM_SENDBUFFER',__LINE__,info)
             GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         !  then allocate the index list: buffer2part
         !----------------------------------------------------------------------
         ldu(1) = (nghostplus - nghost)*ppm_nsublist(topoid)
         CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &           'buffer-to-particles map PPM_BUFFER2PART',__LINE__,info)
             GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         !  flag all ghosts as not yet taken
         !----------------------------------------------------------------------
         lghost(:) = .TRUE.

         !----------------------------------------------------------------------
         !  loop over the subs on the local processor
         !----------------------------------------------------------------------
         DO j=1,ppm_nsublist(topoid)
            !-------------------------------------------------------------------
            !  Get the global ID of the sub
            !-------------------------------------------------------------------
            isub = ppm_isublist(j,topoid) 

            !-------------------------------------------------------------------
            !  Define the extended resize of this sub 
            !-------------------------------------------------------------------
            IF (isymm.GT.0) THEN
               !----------------------------------------------------------------
               !  if we use symmetry ghosts will only be present at the
               !  upper/right part of the sub
               !----------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
               xmini = ppm_min_subs(1,isub,topoid)
               xmaxi = ppm_max_subs(1,isub,topoid) + ghostsize
   
               ymini = ppm_min_subs(2,isub,topoid)
               ymaxi = ppm_max_subs(2,isub,topoid) + ghostsize
   
               IF (ppm_dim.EQ.3) THEN
                  zmini = ppm_min_subs(3,isub,topoid)
                  zmaxi = ppm_max_subs(3,isub,topoid) + ghostsize
               ENDIF 
#else
               xmini = ppm_min_subd(1,isub,topoid)
               xmaxi = ppm_max_subd(1,isub,topoid) + ghostsize
   
               ymini = ppm_min_subd(2,isub,topoid)
               ymaxi = ppm_max_subd(2,isub,topoid) + ghostsize
   
               IF (ppm_dim.EQ.3) THEN
                  zmini = ppm_min_subd(3,isub,topoid)
                  zmaxi = ppm_max_subd(3,isub,topoid) + ghostsize
               ENDIF 
#endif
            ELSE
               !----------------------------------------------------------------
               !  if we do not use symmetry, we have ghost all around
               !----------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
               xmini = ppm_min_subs(1,isub,topoid) - ghostsize
               xmaxi = ppm_max_subs(1,isub,topoid) + ghostsize
   
               ymini = ppm_min_subs(2,isub,topoid) - ghostsize
               ymaxi = ppm_max_subs(2,isub,topoid) + ghostsize
   
               IF (ppm_dim.EQ.3) THEN
                  zmini = ppm_min_subs(3,isub,topoid) - ghostsize
                  zmaxi = ppm_max_subs(3,isub,topoid) + ghostsize
               ENDIF 
#else
               xmini = ppm_min_subd(1,isub,topoid) - ghostsize
               xmaxi = ppm_max_subd(1,isub,topoid) + ghostsize
   
               ymini = ppm_min_subd(2,isub,topoid) - ghostsize
               ymaxi = ppm_max_subd(2,isub,topoid) + ghostsize
   
               IF (ppm_dim.EQ.3) THEN
                  zmini = ppm_min_subd(3,isub,topoid) - ghostsize
                  zmaxi = ppm_max_subd(3,isub,topoid) + ghostsize
               ENDIF 
#endif
            ENDIF 

            !-------------------------------------------------------------------
            !  Reallocate the arrays: ppm_sendbuffers/d and ppm_buffer2part. 
            !  The required size is the current size minus the current use plus 
            !  the maximum increment (nghostplus-nghost).
            !  We perform a grow_preserve to preserve the contents but only to 
            !  reallocate if an increased size is required.
            !-------------------------------------------------------------------
            iopt   = ppm_param_alloc_grow_preserve
            IF (ppm_kind.EQ.ppm_kind_double) THEN
               IF (ASSOCIATED(ppm_sendbufferd)) THEN
                  ldu(1) = SIZE(ppm_sendbufferd) + &
     &                     (nghostplus - nghost - ibuffer)*ppm_dim
               ELSE
                  ldu(1) = (nghostplus - ibuffer)*ppm_dim
               ENDIF 
               CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
               CALL ppm_alloc(ppm_ghost_offsetd,ldu,iopt,info)
            ELSE
               IF (ASSOCIATED(ppm_sendbuffers)) THEN
                  ldu(1) = SIZE(ppm_sendbuffers) + &
     &                     (nghostplus - nghost - ibuffer)*ppm_dim
               ELSE
                  ldu(1) = (nghostplus - nghost - ibuffer)*ppm_dim 
               ENDIF 
               CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
               CALL ppm_alloc(ppm_ghost_offsets,ldu,iopt,info)
            ENDIF
            IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &              'global send buffer PPM_SENDBUFFER',__LINE__,info)
                GOTO 9999
            ENDIF

            IF (ASSOCIATED(ppm_buffer2part)) THEN
               ldu(1) = SIZE(ppm_buffer2part) + nghostplus - nghost - iset
            ELSE
               ldu(1) = nghostplus - nghost - iset
            ENDIF 
            CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
            IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &              'buffer-to-particles map PPM_BUFFER2PART',__LINE__,info)
                GOTO 9999
            ENDIF

            !-------------------------------------------------------------------
            !  loop over the potential ghost particles due to periodicity
            !-------------------------------------------------------------------
            IF (ppm_dim.EQ.2) THEN
               !----------------------------------------------------------------
               !  Two dimensions
               !----------------------------------------------------------------
               DO i=nghost+1,nghostplus
                  !-------------------------------------------------------------
                  !  and check if it is inside the ghost region 
                  !-------------------------------------------------------------
                  IF (xt(1,i).GT.xmini.AND.xt(1,i).LT.xmaxi.AND. &
     &                xt(2,i).GT.ymini.AND.xt(2,i).LT.ymaxi.AND.lghost(i)) THEN
                     !----------------------------------------------------------
                     !  mark the ghost as taken
                     !----------------------------------------------------------
                     lghost(i) = .FALSE.

                     !----------------------------------------------------------
                     !  found one - increment the buffer counter
                     !----------------------------------------------------------
                     iset                  = iset + 1
   
                     !----------------------------------------------------------
                     !  store the ID of the particles
                     !----------------------------------------------------------
                     ppm_buffer2part(iset) = ighost(i)
   
                     !----------------------------------------------------------
                     !  store the particle
                     !----------------------------------------------------------
                     IF (ppm_kind.EQ.ppm_kind_double) THEN
#if    __KIND == __SINGLE_PRECISION
                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd(ibuffer)   = REAL(xt(1,i),        &
     &                       ppm_kind_double)
                         ppm_ghost_offsetd(ibuffer) = REAL(xt_offset(1,i), &
     &                       ppm_kind_double)

                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd(ibuffer)   = REAL(xt(2,i),        &
     &                       ppm_kind_double)
                         ppm_ghost_offsetd(ibuffer) = REAL(xt_offset(2,i), &
     &                       ppm_kind_double)
#else
                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd(ibuffer)   = xt(1,i)
                         ppm_ghost_offsetd(ibuffer) = xt_offset(1,i)

                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd(ibuffer)   = xt(2,i)
                         ppm_ghost_offsetd(ibuffer) = xt_offset(2,i)
#endif
                     ELSE
#if    __KIND == __SINGLE_PRECISION
                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers(ibuffer)   = xt(1,i)
                         ppm_ghost_offsets(ibuffer) = xt_offset(1,i)

                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers(ibuffer)   = xt(2,i)
                         ppm_ghost_offsets(ibuffer) = xt_offset(2,i)
#else
                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers(ibuffer)   = REAL(xt(1,i),        &
     &                       ppm_kind_single)
                         ppm_ghost_offsets(ibuffer) = REAL(xt_offset(1,i), &
     &                       ppm_kind_single)

                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers(ibuffer)   = REAL(xt(2,i),        &
     &                       ppm_kind_single)
                         ppm_ghost_offsets(ibuffer) = REAL(xt_offset(2,i), &
     &                       ppm_kind_single)
#endif
                     ENDIF
                  ENDIF 
               ENDDO
            ELSE
               !----------------------------------------------------------------
               !  Three dimensions
               !----------------------------------------------------------------
               DO i=nghost+1,nghostplus
                  !-------------------------------------------------------------
                  !  and check if it is inside the ghost region 
                  !-------------------------------------------------------------
                  IF (xt(1,i).GT.xmini.AND.xt(1,i).LT.xmaxi.AND. &
     &                xt(2,i).GT.ymini.AND.xt(2,i).LT.ymaxi.AND. & 
     &                xt(3,i).GT.zmini.AND.xt(3,i).LT.zmaxi.AND.lghost(i)) THEN
                     !----------------------------------------------------------
                     !  mark the ghost as taken
                     !----------------------------------------------------------
                     lghost(i) = .FALSE.

                     !----------------------------------------------------------
                     !  found one - increment the buffer counter
                     !----------------------------------------------------------
                     iset                  = iset + 1

                     !----------------------------------------------------------
                     !  store the ID of the particles
                     !----------------------------------------------------------
                     ppm_buffer2part(iset) = ighost(i)

                     !----------------------------------------------------------
                     !  store the particle
                     !----------------------------------------------------------
                     IF (ppm_kind.EQ.ppm_kind_double) THEN
#if    __KIND == __SINGLE_PRECISION
                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd(ibuffer)   = REAL(xt(1,i),        &
     &                      ppm_kind_double)
                         ppm_ghost_offsetd(ibuffer) = REAL(xt_offset(1,i), &
     &                      ppm_kind_double)

                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd(ibuffer)   = REAL(xt(2,i),        &
     &                      ppm_kind_double)
                         ppm_ghost_offsetd(ibuffer) = REAL(xt_offset(2,i), &
     &                      ppm_kind_double)

                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd(ibuffer)   = REAL(xt(3,i),        &
     &                      ppm_kind_double)
                         ppm_ghost_offsetd(ibuffer) = REAL(xt_offset(3,i), &
     &                      ppm_kind_double)
#else
                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd(ibuffer)   = xt(1,i)
                         ppm_ghost_offsetd(ibuffer) = xt_offset(1,i)

                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd(ibuffer)   = xt(2,i)
                         ppm_ghost_offsetd(ibuffer) = xt_offset(2,i)

                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd(ibuffer)   = xt(3,i)
                         ppm_ghost_offsetd(ibuffer) = xt_offset(3,i)
#endif
                     ELSE
#if    __KIND == __SINGLE_PRECISION
                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers(ibuffer)   = xt(1,i)
                         ppm_ghost_offsets(ibuffer) = xt_offset(1,i)

                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers(ibuffer)   = xt(2,i)
                         ppm_ghost_offsets(ibuffer) = xt_offset(2,i)

                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers(ibuffer)   = xt(3,i)
                         ppm_ghost_offsets(ibuffer) = xt_offset(3,i)
#else
                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers(ibuffer)   = REAL(xt(1,i),        &
     &                      ppm_kind_single)
                         ppm_ghost_offsets(ibuffer) = REAL(xt_offset(1,i), &
     &                      ppm_kind_single)

                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers(ibuffer)   = REAL(xt(2,i),        &
     &                      ppm_kind_single)
                         ppm_ghost_offsets(ibuffer) = REAL(xt_offset(2,i), &
     &                      ppm_kind_single)

                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers(ibuffer)   = REAL(xt(3,i),        &
     &                      ppm_kind_single)
                         ppm_ghost_offsets(ibuffer) = REAL(xt_offset(3,i), &
     &                      ppm_kind_single)
#endif
                     ENDIF
                  ENDIF 
               ENDDO
            ENDIF ! 2/3 dimension 
           
         ENDDO ! end of loop over subs on local processor
      ENDIF ! of periodic ghosts

      !-------------------------------------------------------------------------
      !  Update the buffer pointer (ie the current iset ... or the number of
      !  particle we will send to the k-th entry in the icommseq list
      !-------------------------------------------------------------------------
      ppm_psendbuffer(k+1) = iset + 1

      !-------------------------------------------------------------------------
      !  Step 2 and 3:
      !  Notice the ppm_icommseq(1,topoid) points to the local processor, and
      !  we already considered this on step 1, so we now start at k = 2
      !-------------------------------------------------------------------------
      DO k=2,ppm_ncommseq(topoid)
         !----------------------------------------------------------------------
         !  for each processor in the sendrank list we need to send the ghosts
         !  only once, so we flag the lghost as true initially
         !----------------------------------------------------------------------
         lghost(:) = .TRUE.

         !----------------------------------------------------------------------
         !  Get the rank of the processor
         !----------------------------------------------------------------------
         sendrank = ppm_icommseq(k,topoid)
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
         !  only consider positive sendranks
         !----------------------------------------------------------------------
         IF (sendrank.GE.0) THEN
            !-------------------------------------------------------------------
            !  loop over all the subs in the current topology 
            !  this is a bit tedious and perhaps inefficient, but we do not
            !  have the subs belonging to all processors stored for each topoid
            !  sorry ! so we have to perform the calculations each time
            !-------------------------------------------------------------------
            DO j=1,ppm_nsubs(topoid) 
               !----------------------------------------------------------------
               !  consider only those beloging to rank
               !----------------------------------------------------------------
               IF (ppm_subs2proc(j,topoid).EQ.sendrank) THEN
                  !-------------------------------------------------------------
                  !  Define the extended resize of this sub 
                  !-------------------------------------------------------------
                  IF (isymm.GT.0) THEN
                     !----------------------------------------------------------
                     !  if we use symmetry ghosts will only be present at the
                     !  upper/right part of the sub
                     !----------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
                     xmini = ppm_min_subs(1,j,topoid)
                     xmaxi = ppm_max_subs(1,j,topoid) + ghostsize

                     ymini = ppm_min_subs(2,j,topoid)
                     ymaxi = ppm_max_subs(2,j,topoid) + ghostsize

                     IF (ppm_dim.EQ.3) THEN
                        zmini = ppm_min_subs(3,j,topoid)
                        zmaxi = ppm_max_subs(3,j,topoid) + ghostsize
                     ENDIF 
#else
                     xmini = ppm_min_subd(1,j,topoid)
                     xmaxi = ppm_max_subd(1,j,topoid) + ghostsize

                     ymini = ppm_min_subd(2,j,topoid)
                     ymaxi = ppm_max_subd(2,j,topoid) + ghostsize
 
                     IF (ppm_dim.EQ.3) THEN
                        zmini = ppm_min_subd(3,j,topoid)
                        zmaxi = ppm_max_subd(3,j,topoid) + ghostsize
                     ENDIF 
#endif
                  ELSE
                     !----------------------------------------------------------
                     !  if we do not use symmetry, we have ghost all around
                     !----------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
                     xmini = ppm_min_subs(1,j,topoid) - ghostsize
                     xmaxi = ppm_max_subs(1,j,topoid) + ghostsize

                     ymini = ppm_min_subs(2,j,topoid) - ghostsize
                     ymaxi = ppm_max_subs(2,j,topoid) + ghostsize

                     IF (ppm_dim.EQ.3) THEN
                        zmini = ppm_min_subs(3,j,topoid) - ghostsize
                        zmaxi = ppm_max_subs(3,j,topoid) + ghostsize
                     ENDIF 
#else
                     xmini = ppm_min_subd(1,j,topoid) - ghostsize
                     xmaxi = ppm_max_subd(1,j,topoid) + ghostsize

                     ymini = ppm_min_subd(2,j,topoid) - ghostsize
                     ymaxi = ppm_max_subd(2,j,topoid) + ghostsize

                     IF (ppm_dim.EQ.3) THEN
                        zmini = ppm_min_subd(3,j,topoid) - ghostsize
                        zmaxi = ppm_max_subd(3,j,topoid) + ghostsize
                     ENDIF 
#endif
                  ENDIF 

                  !-------------------------------------------------------------
                  !  Reallocate to make sure we have enough memory in the
                  !  ppm_buffer2part and ppm_sendbuffers/d 
                  !-------------------------------------------------------------
                  iopt   = ppm_param_alloc_grow_preserve
                  ldu(1) = ibuffer + nghostplus*ppm_dim
                  IF (ppm_kind.EQ.ppm_kind_double) THEN
                     CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
                     CALL ppm_alloc(ppm_ghost_offsetd,ldu,iopt,info)
                  ELSE
                     CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
                     CALL ppm_alloc(ppm_ghost_offsets,ldu,iopt,info)
                  ENDIF
                  IF (info .NE. 0) THEN
                      info = ppm_error_fatal
                      CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',   &
     &                    'global send buffer PPM_SENDBUFFER',__LINE__,info)
                      GOTO 9999
                  ENDIF

                  ldu(1) = iset + nghostplus
                  CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
                  IF (info .NE. 0) THEN
                      info = ppm_error_fatal
                      CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',   &
     &                   'buffer2particles map PPM_BUFFER2PART',__LINE__,info)
                      GOTO 9999
                  ENDIF

                  !-------------------------------------------------------------
                  !  loop over the potential ghost particles 
                  !-------------------------------------------------------------
                  IF (ppm_dim.EQ.2) THEN
                     !----------------------------------------------------------
                     !  Two dimensions
                     !----------------------------------------------------------
                     DO i=1,nghostplus
                        !-------------------------------------------------------
                        !  and check if it is inside the ghost region 
                        !-------------------------------------------------------
                        IF (xt(1,i).GT.xmini.AND.xt(1,i).LT.xmaxi.AND. &
     &                      xt(2,i).GT.ymini.AND.xt(2,i).LT.ymaxi.AND. &
     &                      lghost(i)) THEN
                           !----------------------------------------------------
                           !  Mark the ghost as taken for this sendrank
                           !----------------------------------------------------
                           lghost(i) = .FALSE.

                           !----------------------------------------------------
                           !  found one - increment the buffer counter
                           !----------------------------------------------------
                           iset                  = iset + 1

                           !----------------------------------------------------
                           !  store the ID of the particles
                           !----------------------------------------------------
                           ppm_buffer2part(iset) = ighost(i)

                           !----------------------------------------------------
                           !  store the particle
                           !----------------------------------------------------
                           IF (ppm_kind.EQ.ppm_kind_double) THEN
#if    __KIND == __SINGLE_PRECISION
                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd(ibuffer)   = REAL(xt(1,i),        &
     &                            ppm_kind_double)
                               ppm_ghost_offsetd(ibuffer) = REAL(xt_offset(1,i), &
     &                            ppm_kind_double)

                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd(ibuffer)   = REAL(xt(2,i),        &
     &                            ppm_kind_double)
                               ppm_ghost_offsetd(ibuffer) = REAL(xt_offset(2,i), &
     &                            ppm_kind_double)
#else
                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd(ibuffer)   = xt(1,i)
                               ppm_ghost_offsetd(ibuffer) = xt_offset(1,i)

                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd(ibuffer)   = xt(2,i)
                               ppm_ghost_offsetd(ibuffer) = xt_offset(2,i)
#endif
                           ELSE
#if    __KIND == __SINGLE_PRECISION
                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers(ibuffer)   = xt(1,i)
                               ppm_ghost_offsets(ibuffer) = xt_offset(1,i)

                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers(ibuffer)   = xt(2,i)
                               ppm_ghost_offsets(ibuffer) = xt_offset(2,i)
#else
                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers(ibuffer)   = REAL(xt(1,i),        &
     &                            ppm_kind_single)
                               ppm_ghost_offsets(ibuffer) = REAL(xt_offset(1,i), &
     &                            ppm_kind_single)

                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers(ibuffer)   = REAL(xt(2,i),   &
     &                            ppm_kind_single)
                               ppm_ghost_offsets(ibuffer) = REAL(xt_offset(2,i), &
     &                            ppm_kind_single)
#endif
                           ENDIF 
                        ENDIF 
                     ENDDO
                  ELSE
                     !----------------------------------------------------------
                     !  Three dimensions
                     !----------------------------------------------------------
                     DO i=1,nghostplus
                        !-------------------------------------------------------
                        !  and check if it is inside the ghost region 
                        !-------------------------------------------------------
                        IF (xt(1,i).GT.xmini.AND.xt(1,i).LT.xmaxi.AND. &
     &                      xt(2,i).GT.ymini.AND.xt(2,i).LT.ymaxi.AND. & 
     &                      xt(3,i).GT.zmini.AND.xt(3,i).LT.zmaxi.AND. &
     &                      lghost(i)) THEN

                           !----------------------------------------------------
                           !  Mark the ghost as taken for this sendrank
                           !----------------------------------------------------
                           lghost(i) = .FALSE.

                           !----------------------------------------------------
                           !  found one - increment the buffer counter
                           !----------------------------------------------------
                           iset                  = iset + 1

                           !----------------------------------------------------
                           !  store the ID of the particles
                           !----------------------------------------------------
                           ppm_buffer2part(iset) = ighost(i)

                           !----------------------------------------------------
                           !  store the particle
                           !----------------------------------------------------
                           IF (ppm_kind.EQ.ppm_kind_double) THEN
#if    __KIND == __SINGLE_PRECISION
                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd(ibuffer)   = REAL(xt(1,i),          &
     &                            ppm_kind_double)
                               ppm_ghost_offsetd(ibuffer) = REAL(xt_offset(1,i), &
     &                            ppm_kind_double)

                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd(ibuffer) = REAL(xt(2,i),          &
     &                            ppm_kind_double)
                               ppm_ghost_offsetd(ibuffer) = REAL(xt_offset(2,i), &
     &                            ppm_kind_double)

                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd(ibuffer)   = REAL(xt(3,i),        &
     &                            ppm_kind_double)
                               ppm_ghost_offsetd(ibuffer) = REAL(xt_offset(3,i), &
     &                            ppm_kind_double)
#else
                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd(ibuffer)   = xt(1,i)
                               ppm_ghost_offsetd(ibuffer) = xt_offset(1,i)
                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd(ibuffer)   = xt(2,i)
                               ppm_ghost_offsetd(ibuffer) = xt_offset(2,i)
                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd(ibuffer)   = xt(3,i)
                               ppm_ghost_offsetd(ibuffer) = xt_offset(3,i)
#endif
                           ELSE
#if    __KIND == __SINGLE_PRECISION
                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers(ibuffer)   = xt(1,i)
                               ppm_ghost_offsets(ibuffer) = xt_offset(1,i)
                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers(ibuffer)   = xt(2,i)
                               ppm_ghost_offsets(ibuffer) = xt_offset(2,i)
                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers(ibuffer)   = xt(3,i)
                               ppm_ghost_offsets(ibuffer) = xt_offset(3,i)
#else
                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers(ibuffer)   = REAL(xt(1,i),   &
     &                            ppm_kind_single)
                               ppm_ghost_offsets(ibuffer) = REAL(xt_offset(1,i), &
     &                            ppm_kind_single)

                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers(ibuffer)   = REAL(xt(2,i),   &
     &                            ppm_kind_single)
                               ppm_ghost_offsets(ibuffer) = REAL(xt_offset(2,i), &
     &                            ppm_kind_single)

                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers(ibuffer)   = REAL(xt(3,i),   &
     &                            ppm_kind_single)
                               ppm_ghost_offsets(ibuffer) = REAL(xt_offset(3,i), &
     &                            ppm_kind_single)
#endif
                           ENDIF 
                        ENDIF 
                     ENDDO
                  ENDIF ! 2/3 dimension 
               ENDIF ! only consider subs belonging to rank
            ENDDO ! loop over all subs in topo
         ENDIF ! skipped negative ranks
         !----------------------------------------------------------------------
         !  Update the buffer pointer (ie the current iset ... or the number of
         !  particle we will send to the k-th entry in the icommseq list
         !----------------------------------------------------------------------
         ppm_psendbuffer(k+1) = iset + 1

      ENDDO ! loop over all processors in commseq

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
         CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ghost_get',     &
     &       'ilist1',__LINE__,info)
      ENDIF

      CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ghost_get',     &
     &       'ilist2',__LINE__,info)
      ENDIF

      CALL ppm_alloc(ighost,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ghost_get',     &
     &       'ighost',__LINE__,info)
      ENDIF

      CALL ppm_alloc(    xt,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ghost_get',     &
     &       'xt',__LINE__,info)
      ENDIF

      CALL ppm_alloc(xt_offset,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ghost_get',     &
     &       'xt',__LINE__,info)
      ENDIF

      CALL ppm_alloc(lghost,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ghost_get',     &
     &       'lghost',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_ghost_get',t0,info)
      RETURN
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_get_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_get_d
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_part_global
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
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
      SUBROUTINE ppm_map_part_global_s(target_topoid,xp,Npart,info,userdef_part2proc)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_global_d(target_topoid,xp,Npart,info,userdef_part2proc)
#endif
      !!! This routine maps the particles onto the given
      !!! topology using a global mapping (i.e. every
      !!! processor communicates with every other).
      !!!
      !!! The first part of the buffer contains the on processor
      !!! data.
      !!!
      !!! The storing of ppm_map_type is not used (yet)
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_map_part_util
      IMPLICIT NONE

#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:)       , INTENT(IN   ) :: xp
      !!! The position of the particles
      INTEGER                        , INTENT(IN   ) :: Npart
      !!! Number of particles (on local processor)
      INTEGER                        , INTENT(IN   ) :: target_topoid
      !!! ID of an existing topology to map onto or `ppm_param_topo_undefined`
      !!! if a non geometric mapping should be used.                           +
      !!! In the latter case ppm will call internally
      !!! `ppm_map_part_equidistributed` to distribute the particles
      !!! proportionaly to processor speeds.
      !!!
      !!! If the user supplies `ppm_param_topo_undefined` an optional parameter
      !!! userdef_part2proc can be used to specify explicitly which particle
      !!! should be moved to which processor.
      INTEGER, DIMENSION(:), OPTIONAL, POINTER       :: userdef_part2proc
      !!! The processor assignment for each particle
      INTEGER                        , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), POINTER :: topo
!       TYPE(ppm_t_topo), POINTER      :: target_topo => NULL()

      REAL(MK) :: t0

      INTEGER, DIMENSION(3) :: ldu
!       INTEGER, DIMENSION(:), POINTER :: bcdef
      INTEGER               :: i,j,k,idom,ipart,nlist1,nlist2
      INTEGER               :: sendrank,recvrank
      INTEGER               :: iopt,iset,ibuffer

      CHARACTER(ppm_char) :: caller='ppm_map_part_global'

      LOGICAL :: valid
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

      ! if there is still some data left in the buffer, warn the user
      IF (ppm_buffer_set .GT. 0) THEN
         fail('Buffer was not empty. Possible loss of data!', &
         & ppm_err_map_incomp,exit_point=no,ppm_error=ppm_error_warning)
      ENDIF

      !-------------------------------------------------------------------------
      !  map to a defined topology
      !-------------------------------------------------------------------------
      IF (target_topoid.NE.ppm_param_topo_undefined) THEN
         topo => ppm_topo(target_topoid)%t

         !-----------------------------------------------------------------------
         !  Save the map type for the subsequent calls (not used yet)
         !-----------------------------------------------------------------------
         ppm_map_type = ppm_param_map_global

         !------------------------------------------------------------------------
         !  Alloc memory for particle lists
         !-----------------------------------------------------------------------
         iopt   = ppm_param_alloc_fit
         ldu(1) = Npart
         CALL ppm_alloc(ilist1,ldu,iopt,info)
         or_fail_alloc("particle list 1 ILIST1",ppm_error=ppm_error_fatal)

         CALL ppm_alloc(ilist2,ldu,iopt,info)
         or_fail_alloc("particle list 2 ILIST2",ppm_error=ppm_error_fatal)

         !-----------------------------------------------------------------------
         !  Allocate memory for the pointer to the buffer; for the global map we
         !  need entries for each processor, thus ldu(1) = ppm_nproc
         !-----------------------------------------------------------------------
         iopt   = ppm_param_alloc_fit
         ldu(1) = ppm_nproc + 1
         CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
         or_fail_alloc("particle send buffer PPM_PSENDBUFFER",ppm_error=ppm_error_fatal)

         iopt   = ppm_param_alloc_fit
         ldu(1) = Npart
         CALL ppm_alloc(part2proc,ldu,iopt,info)
         or_fail_alloc("particles-to-processor map PART2PROC",ppm_error=ppm_error_fatal)

         ldu(1) = Npart
         CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
         or_fail_alloc("buffer-to-particles map PPM_BUFFER2PART",ppm_error=ppm_error_fatal)

         !-----------------------------------------------------------------------
         !  Initialize the particle list
         !-----------------------------------------------------------------------
         FORALL (ipart=1:Npart) ilist1(ipart) = ipart
         nlist1 = Npart

         !-----------------------------------------------------------------------
         !  Then of these exclude the particles that on the current processor
         !-----------------------------------------------------------------------
         IF (ppm_dim.EQ.2) THEN
            !--------------------------------------------------------------------
            !  Loop over the subdomains
            !--------------------------------------------------------------------
  !         DO idom=1,ppm_nsubs(target_topoid)
            DO idom=topo%nsubs,1,-1
               sendrank = topo%sub2proc(idom)

               !-------------------------------------------------------------------
               !  Loop over the remaining particles not yet assigned to a processor
               !-------------------------------------------------------------------
               nlist2 = 0
               DO i=1,nlist1
                  ipart = ilist1(i)
                  !--------------------------------------------------------------
                  !  If the particle is inside the current subdomain, assign it
                  !--------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
                  IF (xp(1,ipart).GE.topo%min_subs(1,idom).AND. &
                  &   xp(1,ipart).LE.topo%max_subs(1,idom).AND. &
                  &   xp(2,ipart).GE.topo%min_subs(2,idom).AND. &
                  &   xp(2,ipart).LE.topo%max_subs(2,idom)) THEN
                     !----------------------------------------------------------
                     !  In the non-periodic case, allow particles that are
                     !  exactly ON an upper EXTERNAL boundary.
                     !----------------------------------------------------------
                     IF ((xp(1,ipart).LT.topo%max_subs(1,idom) .OR. &
                     &  (topo%subs_bc(2,idom).EQ.1           .AND.  &
                     &  topo%bcdef(2).NE.ppm_param_bcdef_periodic)) .AND.  &
                     &  (xp(2,ipart).LT.topo%max_subs(2,idom) .OR.  &
                     &  (topo%subs_bc(4,idom).EQ.1           .AND.  &
                     &  topo%bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
#elif  __KIND == __DOUBLE_PRECISION
                  IF (xp(1,ipart).GE.topo%min_subd(1,idom).AND.   &
                  &   xp(1,ipart).LE.topo%max_subd(1,idom).AND.   &
                  &   xp(2,ipart).GE.topo%min_subd(2,idom).AND.   &
                  &   xp(2,ipart).LE.topo%max_subd(2,idom)) THEN
                     !----------------------------------------------------------
                     !  In the non-periodic case, allow particles that are
                     !  exactly ON an upper EXTERNAL boundary.
                     !----------------------------------------------------------
                     IF ((xp(1,ipart).LT.topo%max_subd(1,idom) .OR. &
                     &  (topo%subs_bc(2,idom).EQ.1           .AND.  &
                     &  topo%bcdef(2).NE.ppm_param_bcdef_periodic)) .AND.  &
                     &  (xp(2,ipart).LT.topo%max_subd(2,idom) .OR.  &
                     &  (topo%subs_bc(4,idom).EQ.1           .AND.  &
                     &  topo%bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
#endif
                        part2proc(ipart) = sendrank
                     ELSE
                        !-------------------------------------------------------
                        !  if not, add it the list of particle to search
                        !-------------------------------------------------------
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

               !-----------------------------------------------------------------
               !  Copy the lists (well, only if nlist2 changed - decreased)
               !-----------------------------------------------------------------
               IF (nlist2.NE.nlist1) THEN
                  nlist1 = nlist2
                  FORALL (i=1:nlist1) ilist1(i)=ilist2(i)
               ENDIF

               !-----------------------------------------------------------------
               !  Exit if the list is empty
               !-----------------------------------------------------------------
               IF (nlist1.EQ.0) EXIT
            ENDDO
         ELSE
            !--------------------------------------------------------------------
            !  Loop over the subdomains (since the first domains are most likely
            !  to be empty, we look backwards to reduce the number of elements in
            !  nlist2 as fast as possible)
            !--------------------------------------------------------------------
  !         DO idom=1,ppm_nsubs(target_topoid)
            DO idom=topo%nsubs,1,-1
               sendrank   = topo%sub2proc(idom)
               !-----------------------------------------------------------------
               !  Loop over the remaining particles
               !-----------------------------------------------------------------
               nlist2 = 0
               DO i=1,nlist1
                  ipart = ilist1(i)
                  !--------------------------------------------------------------
                  !  If the particle is inside the current subdomain, assign it
                  !--------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
                  IF (xp(1,ipart).GE.topo%min_subs(1,idom).AND.   &
                  &   xp(1,ipart).LE.topo%max_subs(1,idom).AND.   &
                  &   xp(2,ipart).GE.topo%min_subs(2,idom).AND.   &
                  &   xp(2,ipart).LE.topo%max_subs(2,idom).AND.   &
                  &   xp(3,ipart).GE.topo%min_subs(3,idom).AND.   &
                  &   xp(3,ipart).LE.topo%max_subs(3,idom)) THEN
                     !----------------------------------------------------------
                     !  In the non-periodic case, allow particles that are
                     !  exactly ON an upper EXTERNAL boundary.
                     !----------------------------------------------------------
                     IF ((xp(1,ipart).LT.topo%max_subs(1,idom) .OR. &
                     &  (topo%subs_bc(2,idom).EQ.1           .AND. &
                     &   topo%bcdef(2).NE. ppm_param_bcdef_periodic))    .AND. &
                     &  (xp(2,ipart).LT.topo%max_subs(2,idom) .OR. &
                     &  (topo%subs_bc(4,idom).EQ.1           .AND. &
                     &   topo%bcdef(4).NE. ppm_param_bcdef_periodic))    .AND. &
                     &  (xp(3,ipart).LT.topo%max_subs(3,idom) .OR. &
                     &  (topo%subs_bc(6,idom).EQ.1           .AND. &
                     &  topo%bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN
#elif  __KIND == __DOUBLE_PRECISION
                  IF (xp(1,ipart).GE.topo%min_subd(1,idom).AND.   &
                  &   xp(1,ipart).LE.topo%max_subd(1,idom).AND.   &
                  &   xp(2,ipart).GE.topo%min_subd(2,idom).AND.   &
                  &   xp(2,ipart).LE.topo%max_subd(2,idom).AND.   &
                  &   xp(3,ipart).GE.topo%min_subd(3,idom).AND.   &
                  &   xp(3,ipart).LE.topo%max_subd(3,idom)) THEN
                     !----------------------------------------------------------
                     !  In the non-periodic case, allow particles that are
                     !  exactly ON an upper EXTERNAL boundary.
                     !----------------------------------------------------------
                     IF ((xp(1,ipart).LT.topo%max_subd(1,idom) .OR. &
                     &  (topo%subs_bc(2,idom).EQ.1           .AND. &
                     &   topo%bcdef(2).NE. ppm_param_bcdef_periodic))    .AND. &
                     &  (xp(2,ipart).LT.topo%max_subd(2,idom) .OR. &
                     &  (topo%subs_bc(4,idom).EQ.1           .AND. &
                     &   topo%bcdef(4).NE. ppm_param_bcdef_periodic))    .AND. &
                     &  (xp(3,ipart).LT.topo%max_subd(3,idom) .OR. &
                     &  (topo%subs_bc(6,idom).EQ.1           .AND. &
                     &   topo%bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN
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
               !-----------------------------------------------------------------
               !  Copy the lists (well, only if nlist2 changed - decreased)
               !-----------------------------------------------------------------
               IF (nlist2.NE.nlist1) THEN
                  nlist1 = nlist2
                  FORALL (i=1:nlist1) ilist1(i)=ilist2(i)
               ENDIF

               !-----------------------------------------------------------------
               !  Exit if the list is empty
               !-----------------------------------------------------------------
               IF (nlist1.EQ.0) EXIT
            ENDDO
         ENDIF

         !-----------------------------------------------------------------------
         !  Here we could check that we sold all the particles, but if we use
         !  Dirichlet BCs we might just want to loose the particles outside of the
         !  domain. So just skip them in the map. This is also consistent with the
         !  map_partial where particles are also lost (there are no periodic sub
         !  images in this case).
         !-----------------------------------------------------------------------
      ELSE
         !-----------------------------------------------------------------------
         ! Non geometric mappings
         !-----------------------------------------------------------------------
         IF (PRESENT(userdef_part2proc)) THEN
            ! just let part2proc point to the user defined mapping
            part2proc => userdef_part2proc
         ENDIF
      ENDIF

      !-----------------------------------------------------------------------
      !  Send the particles to the target processors
      !-----------------------------------------------------------------------
      IF (.NOT.((target_topoid.EQ.ppm_param_topo_undefined).AND. &
      &   (.NOT.PRESENT(userdef_part2proc)))) THEN
        !-----------------------------------------------------------------------
        !  We have a particle processor mapping
        !-----------------------------------------------------------------------


        !-----------------------------------------------------------------------
        !  Store the number of buffer entries (this is the first)
        !-----------------------------------------------------------------------
        ppm_buffer_set = 1

        !-----------------------------------------------------------------------
        !  Allocate memory for the field registers that holds the dimension and
        !  type of the data
        !-----------------------------------------------------------------------
        iopt   = ppm_param_alloc_fit
        ldu(1) = ppm_buffer_set
        CALL ppm_alloc(ppm_buffer_dim ,ldu,iopt,info)
        or_fail_alloc("buffer dimensions PPM_BUFFER_DIM",ppm_error=ppm_error_fatal)

        CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
        or_fail_alloc("buffer types PPM_BUFFER_TYPE",ppm_error=ppm_error_fatal)

        ppm_buffer_dim(ppm_buffer_set)  = ppm_dim
#if    __KIND == __SINGLE_PRECISION
        ppm_buffer_type(ppm_buffer_set) = ppm_kind_single
#else
        ppm_buffer_type(ppm_buffer_set) = ppm_kind_double
#endif

        !-----------------------------------------------------------------------
        !  (Re)allocate memory for the buffer (need not save the contents of the
        !  buffer since this is a global map).
        !-----------------------------------------------------------------------
        iopt   = ppm_param_alloc_fit
        ldu(1) = ppm_dim*Npart
        SELECT CASE (ppm_kind)
        CASE (ppm_kind_double)
           CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
        CASE DEFAULT
           CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
        END SELECT
        or_fail_alloc("global send buffer PPM_SENDBUFFER",ppm_error=ppm_error_fatal)

        !-----------------------------------------------------------------------
        !  Allocate memory for the sendlist
        !-----------------------------------------------------------------------
        ppm_nsendlist = ppm_nproc
        ppm_nrecvlist = ppm_nproc
        ldu(1)        = ppm_nsendlist
        CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
        or_fail_alloc("send list PPM_ISENDLIST",ppm_error=ppm_error_fatal)

        CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
        or_fail_alloc("receive list PPM_IRECVLIST",ppm_error=ppm_error_fatal)

        !-----------------------------------------------------------------------
        !  Initialize the particle lists
        !-----------------------------------------------------------------------
        nlist1 = Npart
        FORALL (ipart=1:Npart) ilist1(ipart) = ipart

        !-----------------------------------------------------------------------
        !  loop over all processors, starting with the processor itself
        !-----------------------------------------------------------------------
        sendrank           = ppm_rank - 1
        recvrank           = ppm_rank + 1
        ppm_psendbuffer(1) = 1
        ppm_nsendlist      = 0
        ppm_nrecvlist      = 0
        iset               = 0
        ibuffer            = 0

        DO i=1,ppm_nproc
           !--------------------------------------------------------------------
           !  compute the next processor
           !--------------------------------------------------------------------
           sendrank = sendrank + 1
           IF (sendrank.GT.ppm_nproc-1) sendrank = sendrank - ppm_nproc
           recvrank = recvrank - 1
           IF (recvrank.LT.          0) recvrank = recvrank + ppm_nproc

           !--------------------------------------------------------------------
           !  Store the processor to which we will send to
           !--------------------------------------------------------------------
           ppm_nsendlist                = ppm_nsendlist + 1
           ppm_isendlist(ppm_nsendlist) = sendrank

           !--------------------------------------------------------------------
           !  Store the processor to which we will recv from
           !--------------------------------------------------------------------
           ppm_nrecvlist                = ppm_nrecvlist + 1
           ppm_irecvlist(ppm_nrecvlist) = recvrank

           !--------------------------------------------------------------------
           !  Initialize the buffer count
           !--------------------------------------------------------------------
           nlist2 = 0

           IF (ppm_dim .EQ. 2) THEN
              DO j=1,nlist1
                 ipart = ilist1(j)
                 IF (part2proc(ipart).EQ.sendrank) THEN
                    !-----------------------------------------------------------
                    !  increment the buffer counter
                    !-----------------------------------------------------------
                    iset = iset + 1

                    !-----------------------------------------------------------
                    !  Store the id of the particle
                    !-----------------------------------------------------------
                    ppm_buffer2part(iset) = ipart

                    !-----------------------------------------------------------
                    !  Store the particle
                    !-----------------------------------------------------------
                    IF (ppm_kind .EQ. ppm_kind_double) THEN
#if    __KIND == __SINGLE_PRECISION
                       ibuffer                  = ibuffer + 1
                       ppm_sendbufferd(ibuffer) = REAL(xp(1,ipart),ppm_kind_double)
                       ibuffer                  = ibuffer + 1
                       ppm_sendbufferd(ibuffer) = REAL(xp(2,ipart),ppm_kind_double)
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
                       ppm_sendbuffers(ibuffer) = REAL(xp(1,ipart),ppm_kind_single)
                       ibuffer                  = ibuffer + 1
                       ppm_sendbuffers(ibuffer) = REAL(xp(2,ipart),ppm_kind_single)
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
                    !-----------------------------------------------------------
                    !  increment the buffer counter
                    !-----------------------------------------------------------
                    iset = iset + 1

                    !-----------------------------------------------------------
                    !  Store the id of the particle
                    !-----------------------------------------------------------
                    ppm_buffer2part(iset) = ipart

                    !-----------------------------------------------------------
                    !  Store the particle
                    !-----------------------------------------------------------
                    IF (ppm_kind .EQ. ppm_kind_double) THEN
#if    __KIND == __SINGLE_PRECISION
                       ibuffer                  = ibuffer + 1
                       ppm_sendbufferd(ibuffer) = REAL(xp(1,ipart),ppm_kind_double)
                       ibuffer                  = ibuffer + 1
                       ppm_sendbufferd(ibuffer) = REAL(xp(2,ipart),ppm_kind_double)
                       ibuffer                  = ibuffer + 1
                       ppm_sendbufferd(ibuffer) = REAL(xp(3,ipart),ppm_kind_double)
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
                       ppm_sendbuffers(ibuffer) = REAL(xp(1,ipart),ppm_kind_single)
                       ibuffer                  = ibuffer + 1
                       ppm_sendbuffers(ibuffer) = REAL(xp(2,ipart),ppm_kind_single)
                       ibuffer                  = ibuffer + 1
                       ppm_sendbuffers(ibuffer) = REAL(xp(3,ipart),ppm_kind_single)
#endif
                    ENDIF
                 ELSE
                    nlist2         = nlist2 + 1
                    ilist2(nlist2) = ipart
                 ENDIF
              ENDDO
           ENDIF

           !--------------------------------------------------------------------
           !  Update the buffer pointer
           !--------------------------------------------------------------------
           ppm_psendbuffer(i+1) = iset + 1

           !--------------------------------------------------------------------
           !  Swap the lists
           !--------------------------------------------------------------------
           nlist1 = nlist2
           DO j=1,nlist1
              ilist1(j) = ilist2(j)
           ENDDO
        ENDDO

        !-----------------------------------------------------------------------
        !  Everything has to go (should not be necessary!)
        !  the line above:
        !
        !      ppm_psendbuffer(i+1) = iset + 1
        !
        !  should take care of this
        !-----------------------------------------------------------------------
        ! PC: killed this for unbounded domains where the particles can leave the
        ! the domain
        ! ppm_psendbuffer(ppm_nproc+1) = npart + 1

        !-----------------------------------------------------------------------
        !  Store the current size of the buffer
        !-----------------------------------------------------------------------
        ppm_nsendbuffer = ibuffer

        !-----------------------------------------------------------------------
        !  Deallocate the memory for the lists
        !-----------------------------------------------------------------------
        iopt = ppm_param_dealloc
        CALL ppm_alloc(ilist1,ldu,iopt,info)
        or_fail_dealloc('particle list 1 ILIST1',exit_point=no)

        CALL ppm_alloc(ilist2,ldu,iopt,info)
        or_fail_dealloc('particle list 2 ILIST2',exit_point=no)

        IF (.NOT. PRESENT(userdef_part2proc)) THEN
           CALL ppm_alloc(part2proc,ldu,iopt,info)
           or_fail_dealloc('particle-to-processor map PART2PROC',exit_point=no)
        ENDIF
      ELSE
        !-------------------------------------------------------------------
        ! equidistributed mapping, handled by ppm_map_part_eqdistrib
        !-------------------------------------------------------------------
        CALL ppm_map_part_eqdistrib(xp,Npart,info)
        or_fail('Failed to equidistribute particles',ppm_error=ppm_error_fatal)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (Npart .LT. 0) THEN
           fail('Npart must be >0',exit_point=8888)
        ENDIF
        IF (target_topoid .NE. ppm_param_topo_undefined) THEN
            CALL ppm_check_topoid(target_topoid,valid,info)
            IF (.NOT. valid) THEN
               fail('target_target_topoid is not valid',exit_point=8888)
            ENDIF
        ENDIF
        IF (PRESENT(userdef_part2proc)) THEN
            IF (.NOT. ASSOCIATED(userdef_part2proc)) THEN
               fail('userdef_part2proc has to be associated',exit_point=8888)
            ENDIF
            IF (SIZE(userdef_part2proc) .LT. Npart) THEN
               fail('userdef_part2proc has to be at least Npart large',exit_point=8888)
            ENDIF
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_global_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_global_d
#endif

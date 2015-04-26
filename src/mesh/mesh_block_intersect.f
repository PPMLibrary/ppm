      SUBROUTINE equi_mesh_block_intersect(this,to_mesh,isub,jsub,offsetb, &
      &          nsendlistb,isendfromsubb,isendtosubb,isendpatchidb,       &
      &          isendblkstartb,isendblksizeb,ioffsetb,info,lsymm,ghostsize)
      !!! This routine determines, for each patch, common mesh blocks
      !!! (intersections) of two subpatches, on different subdomains,
      !!! possibly shifted and/or extended with ghostlayers.
      !!! Start and size of all blocks found are stored and returned. This
      !!! routine is used by `equi_mesh_map_global` and
      !!! `equi_mesh_map_ghost_init`.
      !!!
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_topo_typedef, ONLY : ppm_topo
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_t_equi_mesh)                         :: this
      !!! Source mesh
      CLASS(ppm_t_equi_mesh_),         INTENT(IN   ) :: to_mesh
      !!! Target mesh
      INTEGER,                         INTENT(IN   ) :: isub
      !!! Source sub (from which data will be sent) in global numbering
      INTEGER,                         INTENT(IN   ) :: jsub
      !!! Destination sub in global numbering
      INTEGER, DIMENSION(:),           INTENT(IN   ) :: offsetb
      !!! Shift offset for the source subdomain. Index: 1:ppm_dim
      INTEGER,                         INTENT(INOUT) :: nsendlistb
      !!! Number of mesh blocks in the lists so far
      INTEGER, DIMENSION(:),           POINTER       :: isendfromsubb
      !!! Global sub index of sources
      INTEGER, DIMENSION(:),           POINTER       :: isendtosubb
      !!! Global sub index of targets
      INTEGER, DIMENSION(:,:),         POINTER       :: isendpatchidb
      !!! Global patch index of data patches
      INTEGER, DIMENSION(:,:),         POINTER       :: isendblkstartb
      !!! Start of the mesh blocks in the global mesh. 1st index:
      !!! 1:ppm_dim, 2nd: meshblock.
      INTEGER, DIMENSION(:,:),         POINTER       :: isendblksizeb
      !!! Size of the meshblocks in numbers of grid points. 1st index:
      !!! 1:ppm_dim, 2nd: meshblock
      INTEGER, DIMENSION(:,:),         POINTER       :: ioffsetb
      !!! Meshblock offset for periodic images.
      !!!
      !!! 1st index: 1:ppm_dim                                                 +
      !!! 2nd index: meshblock
      !!!
      !!! To be added to isendblkstart in order to get irecvblkstart on the
      !!! target subdomain.
      INTEGER,                         INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      LOGICAL, DIMENSION(3), OPTIONAL, INTENT(IN   ) :: lsymm
      !!! Use symmetry and chop last point
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: ghostsize
      !!! Size of the ghost layer in numbers of grid points in all space
      !!! dimensions (1...ppm_dim).

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(2)       :: ldu
      INTEGER, DIMENSION(ppm_dim) :: iblockstart,nblocksize
      INTEGER                     :: iblockstopk,k,iopt,isize,ipatch
      INTEGER                     :: fromistartk,toistartk
      INTEGER                     :: fromiendk,toiendk
      INTEGER, DIMENSION(ppm_dim) :: chop
      INTEGER, DIMENSION(ppm_dim) :: ghostsize_

      LOGICAL :: dosend

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      start_subroutine("mesh_block_intersect")

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine whether to send last point too (use symmetry)
      !-------------------------------------------------------------------------
      IF (PRESENT(lsymm)) THEN
         chop=MERGE(1,0,lsymm)
      ELSE
         chop=0
      END IF

      IF (PRESENT(ghostsize)) THEN
         ghostsize_=MIN(ghostsize,to_mesh%ghostsize)
      ELSE
         ghostsize_=to_mesh%ghostsize
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine current minimal common length of lists
      !-------------------------------------------------------------------------
      isize = MIN(SIZE(isendfromsubb,1),  &
      &           SIZE(isendtosubb,1),    &
      &           SIZE(isendpatchidb,2),  &
      &           SIZE(isendblkstartb,2), &
      &           SIZE(isendblksizeb,2),  &
      &           SIZE(ioffsetb,2))

      !-------------------------------------------------------------------------
      !  Determine intersecting mesh blocks
      !-------------------------------------------------------------------------
      ! loop over each subpatch of the target sub
      IF (ppm_debug.GT.2) THEN
         stdout("FOR source_sub = ",isub," target_sub",jsub)
         ! yaser: subpach_by_sub should contain global index
         ! so I would use jsub instead of isub_loc
         stdout("   nb subpatches = ",'to_mesh%subpatch_by_sub(jsub)%nsubpatch')
      ENDIF
      ! yaser: subpach_by_sub should contain global index
      ! so I would use jsub instead of isub_loc
      !print*,ppm_rank,isub,jsub,to_mesh%subpatch_by_sub(jsub)%nsubpatch,to_mesh%subpatch_by_sub(isub)%nsubpatch
      DO ipatch=1,to_mesh%subpatch_by_sub(jsub)%nsubpatch
         dosend = .TRUE.
         SELECT TYPE(p => to_mesh%subpatch_by_sub(jsub)%vec(ipatch)%t)
         TYPE IS (ppm_t_subpatch)
            !intersect the patch that this subpatch belongs to
            ! (coordinates istart_p and iend_p) with the extended target
            ! sub and with the source sub
            IF (ppm_debug.GT.2) THEN
               stdout("---------------------------------")
               stdout("INTERSECTING: source_sub = ",isub," target_sub",jsub)
               stdout("  OFFSET = ",offsetb)
               stdout("  CHOP = ",chop)
               stdout(" patch : ")
               stdout("     ",'p%istart_p')
               stdout("     ",'p%iend_p')
               stdout(" subpatch : ")
               stdout("     ",'p%istart')
               stdout("     ",'p%iend')
               stdout(" ghostsizes : ")
               stdout("     ",'p%ghostsize(1)','p%ghostsize(3)')
               stdout("     ",'p%ghostsize(2)','p%ghostsize(4)')
               stdout(" source sub : ")
               stdout("     ",'this%istart(1:ppm_dim,isub)')
               stdout("     ",'this%iend(1:ppm_dim,isub)')
               stdout(" offset source sub : ")
               stdout("     ",'this%istart(1:ppm_dim,isub)+offsetb(1:ppm_dim)')
               stdout("     ",'this%iend(1:ppm_dim,isub)+offsetb(1:ppm_dim)-chop(1:ppm_dim)')
               stdout(" target sub : ")
               stdout("     ",'to_mesh%istart(1:ppm_dim,jsub)')
               stdout("     ",'to_mesh%iend(1:ppm_dim,jsub)')
               stdout(" extended target sub : ")
               stdout("     ",'to_mesh%istart(1,jsub)-this%ghostsize(1)',&
                              'to_mesh%istart(2,jsub)-this%ghostsize(1)')
               IF (ppm_dim.GT.2) THEN
                  stdout("     ",'to_mesh%istart(3,jsub)-this%ghostsize(3)')
               END IF
               stdout("     ",'to_mesh%iend(1,jsub)+this%ghostsize(2)',&
                              'to_mesh%iend(2,jsub)+this%ghostsize(2)')
               IF (ppm_dim.GT.2) THEN
                  stdout("     ",'to_mesh%iend(3,jsub)+this%ghostsize(3)')
               END IF
            ENDIF

            !No intersection if the patch is not on the target subdomain
            DO k=1,ppm_dim
               IF (p%istart_p(k).GT.to_mesh%iend(k,jsub)) THEN
                  dosend = .FALSE.
                  EXIT
               ENDIF
               IF (p%iend_p(k).LT.to_mesh%istart(k,jsub)) THEN
                  dosend = .FALSE.
                  EXIT
               ENDIF

               fromistartk = this%istart(k,isub)+offsetb(k)
               fromiendk   = this%iend(k,isub)  +offsetb(k)-chop(k)

               toistartk = p%istart_p(k)
               toiendk   = p%iend_p(k)

               iblockstart(k) = MAX(fromistartk,toistartk)
               iblockstopk    = MIN(fromiendk,toiendk)

               ! (Here we would need to use the ghostlayer sizes of the target
               ! subpatch, which lives on another subdomain and which we
               ! therefore cannot access as it could be on another proc.
               ! This information could be stored, but instead we just
               ! use the global ghostsize and reject the cases that fail.
               ! The latter are those for which the source patch does not
               ! have real mesh nodes on the target sub).
               toistartk = to_mesh%istart(k,jsub)-ghostsize_(k)
               toiendk   = to_mesh%iend(k,jsub)  +ghostsize_(k)

               iblockstart(k) = MAX(iblockstart(k),toistartk)
               iblockstopk    = MIN(iblockstopk,   toiendk)

               !size of the block (+1 if mesh nodes at the boundary
               ! between blocks are NOT duplicated, and should thus be
               ! counted as ghosts)
               nblocksize(k)  = iblockstopk-iblockstart(k)+1
               ! do not send if there is nothing to be sent
               IF (nblocksize(k).LT.1) THEN
                  dosend = .FALSE.
                  EXIT
               ENDIF
            ENDDO
            !Would that make a difference in speed?
            !dosend = ALL(nblocksize.GE.1)

            IF (ppm_debug.GT.2) THEN
               stdout("iblockstart = ",iblockstart)
               stdout("nblocksize = ",nblocksize)
               stdout("dosend = ",dosend)
               stdout("---------------------------------")
            ENDIF

            IF (dosend) THEN
               nsendlistb = nsendlistb + 1
               IF (nsendlistb .GT. isize) THEN
                  !-----------------------------------------------------------
                  !  Grow memory for the sendlists
                  !-----------------------------------------------------------
                  iopt   = ppm_param_alloc_grow_preserve
                  isize  = MAX(2*isize,nsendlistb)
                  ldu(1) = isize
                  CALL ppm_alloc(isendfromsubb,ldu,iopt,info)
                  or_fail_alloc("isendfromsub")

                  CALL ppm_alloc(isendtosubb,ldu,iopt,info)
                  or_fail_alloc("isendtosub")

                  ldu(1) = ppm_dim
                  ldu(2) = isize
                  CALL ppm_alloc(isendpatchidb,ldu,iopt,info)
                  or_fail_alloc("isendpatchid")

                  CALL ppm_alloc(isendblkstartb,ldu,iopt,info)
                  or_fail_alloc("isendblkstart")

                  CALL ppm_alloc(isendblksizeb,ldu,iopt,info)
                  or_fail_alloc("isendblksize")

                  CALL ppm_alloc(ioffsetb,ldu,iopt,info)
                  or_fail_alloc("ioffset")
               ENDIF
               !---------------------------------------------------------------
               !  send block (isendblkstart...isendblkstart+nblocksize-1)
               !  from sub isub to sub jsub
               !---------------------------------------------------------------
               ! global sub index of local sub where the blocks come from
               isendfromsubb(nsendlistb)           = isub
               ! global sub index of other sub where the blocks go to
               isendtosubb(nsendlistb)             = jsub
               ! global patch index of which the blocks belong
               isendpatchidb(1:ppm_dim,nsendlistb) = p%istart_p(1:ppm_dim)
               ! start of mesh block in global mesh
               isendblkstartb(1:ppm_dim,nsendlistb)= iblockstart(1:ppm_dim)-offsetb(1:ppm_dim)
               ! size of mesh block
               isendblksizeb(1:ppm_dim,nsendlistb) = nblocksize(1:ppm_dim)
               ! mesh block offset for periodic images
               ioffsetb(1:ppm_dim,nsendlistb)      = offsetb(1:ppm_dim)
            ENDIF
         END SELECT
      ENDDO

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      end_subroutine()

      RETURN
      CONTAINS
      SUBROUTINE check
          IF (SIZE(offsetb,1) .LT. ppm_dim) THEN
             fail('offset must be at least of length ppm_dim',exit_point=8888)
          ENDIF
          IF ((isub.LE.0).OR.(isub.GT.ppm_topo(this%topoid)%t%nsubs)) THEN
             fail('isub out of range',exit_point=8888)
          ENDIF
          IF ((jsub.LE.0).OR.(jsub.GT.ppm_topo(to_mesh%topoid)%t%nsubs)) THEN
             fail('jsub out of range',exit_point=8888)
          ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE equi_mesh_block_intersect

      SUBROUTINE equi_mesh_block_intersect(this,to_mesh,isub,jsub,offset, &
      &          nsendlist,isendfromsub,isendtosub,isendpatchid,          &
      &          isendblkstart,isendblksize,ioffset,info,lsymm)
      !!! This routine determines, for each patch, common mesh blocks
      !!! (intersections) of two subpatches, on different subdomains,
      !!! possibly shifted and/or extended with ghostlayers.
      !!! Start and size of all blocks found are stored and returned. This
      !!! routine is used by `ppm_map_field_global` and
      !!! `ppm_map_field_ghost_init`.
      !!!
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data_mesh
      USE ppm_module_topo_typedef
      USE ppm_module_check_id
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_t_equi_mesh),   INTENT(IN   ) :: this
      !!! Source mesh
      CLASS(ppm_t_equi_mesh_),  INTENT(IN   ) :: to_mesh
      !!! Target mesh
      INTEGER,                  INTENT(IN   ) :: isub
      !!! Source sub (from which data will be sent) in global numbering
      INTEGER,                  INTENT(IN   ) :: jsub
      !!! Destination sub in global numbering
      INTEGER, DIMENSION(:),    INTENT(IN   ) :: offset
      !!! Shift offset for the source subdomain. Index: 1:ppm_dim
      INTEGER, DIMENSION(:),    POINTER       :: isendfromsub
      !!! Global sub index of sources
      INTEGER, DIMENSION(:),    POINTER       :: isendtosub
      !!! Global sub index of targets
      INTEGER, DIMENSION(:,:),  POINTER       :: isendpatchid
      !!! Global patch index of data patches
      INTEGER, DIMENSION(:,:),  POINTER       :: ioffset
      !!! Meshblock offset for periodic images.
      !!!
      !!! 1st index: 1:ppm_dim                                                 +
      !!! 2nd index: meshblock
      !!!
      !!! To be added to isendblkstart in order to get irecvblkstart on the
      !!! target subdomain.
      INTEGER, DIMENSION(:,:),  POINTER       :: isendblkstart
      !!! Start of the mesh blocks in the global mesh. 1st index:
      !!! 1:ppm_dim, 2nd: meshblock.
      INTEGER, DIMENSION(:,:),  POINTER       :: isendblksize
      !!! Size of the meshblocks in numbers of grid points. 1st index:
      !!! 1:ppm_dim, 2nd: meshblock
      INTEGER,                  INTENT(INOUT) :: nsendlist
      !!! Number of mesh blocks in the lists so far
      INTEGER,                  INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      LOGICAL, DIMENSION(3),    INTENT(IN   ), OPTIONAL :: lsymm
      !!! Use symmetry and chop last point
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(ppm_t_subpatch),  POINTER :: p => NULL()

      INTEGER, DIMENSION(2)       :: ldu
      INTEGER, DIMENSION(ppm_dim) :: iblockstart,nblocksize
      INTEGER                     :: iblockstopk,k,iopt,pdim,isize,ipatch
      INTEGER                     :: fromistartk,toistartk
      INTEGER                     :: fromiendk,toiendk
      INTEGER, DIMENSION(ppm_dim) :: chop
      LOGICAL                     :: dosend
      LOGICAL                     :: valid

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      start_subroutine("mesh_block_intersect")

      pdim = ppm_dim

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

      !-------------------------------------------------------------------------
      !  Determine current minimal common length of lists
      !-------------------------------------------------------------------------
      isize = MIN(SIZE(isendfromsub,1), SIZE(isendtosub,1),  SIZE(isendpatchid,2),&
      &           SIZE(isendblkstart,2),SIZE(isendblksize,2),SIZE(ioffset,2))

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
                 stdout("  OFFSET = ",offset)
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
                 stdout("     ",'this%istart(1:pdim,isub)')
                 stdout("     ",'this%iend(1:pdim,isub)')
                 stdout(" offset source sub : ")
                 stdout("     ",'this%istart(1:pdim,isub)+offset(1:pdim)')
                 stdout("     ",'this%iend(1:pdim,isub)+offset(1:pdim)-chop(1:pdim)')
                 stdout(" target sub : ")
                 stdout("     ",'to_mesh%istart(1:pdim,jsub)')
                 stdout("     ",'to_mesh%iend(1:pdim,jsub)')
                 stdout(" extended target sub : ")
                 stdout("     ",'to_mesh%istart(1,jsub)-this%ghostsize(1)',&
                                'to_mesh%istart(2,jsub)-this%ghostsize(1)')
                 IF (pdim.GT.2) THEN
                    stdout("     ",'to_mesh%istart(3,jsub)-this%ghostsize(3)')
                 END IF
                 stdout("     ",'to_mesh%iend(1,jsub)+this%ghostsize(2)',&
                                'to_mesh%iend(2,jsub)+this%ghostsize(2)')
                 IF (pdim.GT.2) THEN
                    stdout("     ",'to_mesh%iend(3,jsub)+this%ghostsize(3)')
                 END IF
              ENDIF

              !No intersection if the patch is not on the target subdomain
              DO k=1,pdim
                 fromistartk = this%istart(k,isub)+offset(k)
                 fromiendk   = this%iend(k,isub)  +offset(k)-chop(k)

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
                 toistartk = to_mesh%istart(k,jsub)-to_mesh%ghostsize(k)
                 toiendk   = to_mesh%iend(k,jsub)  +to_mesh%ghostsize(k)

                 iblockstart(k) = MAX(iblockstart(k),toistartk)
                 iblockstopk    = MIN(iblockstopk,   toiendk)

                 !size of the block (+1 if mesh nodes at the boundary
                 ! between blocks are NOT duplicated, and should thus be
                 ! counted as ghosts)
                 nblocksize(k)  = iblockstopk-iblockstart(k)+1
                 ! do not send if there is nothing to be sent
                 IF (nblocksize(k).LT.1) dosend = .FALSE.

                 IF (p%istart_p(k).GT.to_mesh%iend(k,jsub)) dosend = .FALSE.
                 IF (p%iend_p(k).LT.to_mesh%istart(k,jsub)) dosend = .FALSE.
              ENDDO
              !Would that make a difference in speed?
              !dosend = ALL(nblocksize.GE.1)

              IF (ppm_debug.GT.2) THEN
                 stdout("iblockstart = ",iblockstart)
                 stdout("nblockssize = ",nblocksize)
                 stdout("dosend = ",dosend)
                 stdout("---------------------------------")
              ENDIF

              IF (dosend) THEN
                 nsendlist = nsendlist + 1
                 IF (nsendlist .GT. isize) THEN
                    !-----------------------------------------------------------
                    !  Grow memory for the sendlists
                    !-----------------------------------------------------------
                    isize  = MAX(2*isize,nsendlist)
                    iopt   = ppm_param_alloc_grow_preserve
                    ldu(1) = isize
                    CALL ppm_alloc(isendfromsub,ldu,iopt,info)
                    or_fail_alloc("isendfromsub")

                    CALL ppm_alloc(isendtosub,ldu,iopt,info)
                    or_fail_alloc("isendtosub")

                    ldu(1) = pdim
                    ldu(2) = isize
                    CALL ppm_alloc(isendpatchid,ldu,iopt,info)
                    or_fail_alloc("isendpatchid")

                    CALL ppm_alloc(isendblkstart,ldu,iopt,info)
                    or_fail_alloc("isendblkstart")

                    CALL ppm_alloc(isendblksize,ldu,iopt,info)
                    or_fail_alloc("isendblksize")

                    CALL ppm_alloc(ioffset,ldu,iopt,info)
                    or_fail_alloc("ioffset")
                 ENDIF
                 !---------------------------------------------------------------
                 !  send block (isendblkstart...isendblkstart+nblocksize-1)
                 !  from sub isub to sub jsub
                 !---------------------------------------------------------------
                 ! global sub index of local sub where the blocks come from
                 isendfromsub(nsendlist)        = isub
                 ! global sub index of other sub where the blocks go to
                 isendtosub(nsendlist)          = jsub
                 ! global patch index of which the blocks belong
                 isendpatchid(1:pdim,nsendlist) = p%istart_p(1:pdim)
                 ! start of mesh block in global mesh
                 isendblkstart(1:pdim,nsendlist)= iblockstart(1:pdim) - &
                 &                                offset(1:pdim)
                 ! size of mesh block
                 isendblksize(1:pdim,nsendlist) = nblocksize(1:pdim)
                 ! mesh block offset for periodic images
                 ioffset(1:pdim,nsendlist)      = offset(1:pdim)
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
          IF (SIZE(offset,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
              & 'offset must be at least of length ppm_dim',__LINE__,info)
              GOTO 8888
          ENDIF
          IF ((isub.LE.0).OR.(isub.GT.ppm_topo(this%topoid)%t%nsubs)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
              & 'isub out of range',__LINE__,info)
              GOTO 8888
          ENDIF
          IF ((jsub.LE.0).OR.(jsub.GT.ppm_topo(to_mesh%topoid)%t%nsubs)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
              & 'jsub out of range',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE equi_mesh_block_intersect

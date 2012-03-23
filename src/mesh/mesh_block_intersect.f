      SUBROUTINE equi_mesh_block_intersect(this,to_mesh,isub,jsub,offset,&
              ghostsize,nsendlist,isendfromsub,isendtosub,&
              isendblkstart,isendblksize,ioffset,info,lsymm)
      !!! This routine determines common mesh blocks (intersections) of
      !!! two subs, possibly shifted and/or extended with ghostlayers.
      !!! Start and size of all blocks found are stored and returned. This
      !!! routine is used by `ppm_map_field_global` and
      !!! `ppm_map_field_ghost_init`.
      !!!
      !!! [WARNING]
      !!! This does not check whether the from_mesh and the to_mesh are
      !!! compatible (i.e. have the same Nm). This check must be done by
      !!! the calling routine. Furthermore, the 5 output lists need to be
      !!! allocated by the calling routine. They are just grown here if needed.
      !!! The reason is that this routine is typically called from inside a
      !!! loop and we do not want to redo the checks at every iteration.
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
      CLASS(ppm_t_equi_mesh)                  :: this
      !!! Source mesh
      CLASS(ppm_t_equi_mesh_)                 :: to_mesh
      !!! Target mesh
      INTEGER                 , INTENT(IN   ) :: isub
      !!! Source sub (from which data will be sent) in global numbering
      INTEGER                 , INTENT(IN   ) :: jsub
      !!! Destination sub in global numbering
      INTEGER, DIMENSION(:)   , INTENT(IN   ) :: ghostsize
      !!! Size of the ghost layer in numbers of grid points in all space
      !!! dimensions (1...ppm_dim). The target subdomain will be enlarged
      !!! by this.
      INTEGER, DIMENSION(:)   , INTENT(IN   ) :: offset
      !!! Shift offset for the source subdomain. Index: 1:ppm_dim
      INTEGER, DIMENSION(:)   , POINTER       :: isendfromsub
      !!! Global sub index of sources
      INTEGER, DIMENSION(:)   , POINTER       :: isendtosub
      !!! Global sub index of targets
      INTEGER, DIMENSION(:,:) , POINTER       :: ioffset
      !!! Meshblock offset for periodic images.
      !!!
      !!! 1st index: 1:ppm_dim                                                 +
      !!! 2nd index: meshblock
      !!!
      !!! To be added to isendblkstart in order to get irecvblkstart on the
      !!! target subdomain.
      INTEGER, DIMENSION(:,:) , POINTER       :: isendblkstart
      !!! Start of the mesh blocks in the global mesh. 1st index:
      !!! 1:ppm_dim, 2nd: meshblock.
      INTEGER, DIMENSION(:,:) , POINTER       :: isendblksize
      !!! Size of the meshblocks in numbers of grid points. 1st index:
      !!! 1:ppm_dim, 2nd: meshblock
      INTEGER                 , INTENT(INOUT) :: nsendlist
      !!! Number of mesh blocks in the lists so far
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      LOGICAL, DIMENSION(3)   , INTENT(IN   ), OPTIONAL :: lsymm
      !!! Use symmetry and chop last point
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)            :: ldu
      INTEGER, DIMENSION(ppm_dim)      :: iblockstart,nblocksize
      INTEGER                          :: iblockstopk,k,iopt,pdim,isize
      INTEGER                          :: fromistartk,toistartk
      INTEGER                          :: fromnnodesk,tonnodesk
      INTEGER, DIMENSION(3)            :: chop
      LOGICAL                          :: dosend
      LOGICAL                          :: valid
      TYPE(ppm_t_topo),      POINTER   :: from_topo => NULL()
      TYPE(ppm_t_topo),      POINTER   :: to_topo   => NULL()
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

      from_topo => ppm_topo(this%topoID)%t
      to_topo => ppm_topo(to_mesh%topoID)%t

      !-------------------------------------------------------------------------
      !  Determine whether to send last point too (use symmetry)
      !-------------------------------------------------------------------------
      chop = 0
      IF (PRESENT(lsymm)) THEN
         IF (lsymm(1).EQV..TRUE.) THEN
            chop(1) = 1
         END IF
         IF (lsymm(2).EQV..TRUE.) THEN
            chop(2) = 1
         END IF
         IF (pdim .GT. 2) THEN
            IF (lsymm(3).EQV..TRUE.) THEN
               chop(3) = 1
            END IF
         END IF
      END IF

      !-------------------------------------------------------------------------
      !  Determine current minimal common length of lists
      !-------------------------------------------------------------------------
      isize = MIN(SIZE(isendfromsub,1),SIZE(isendtosub,1),        &
     &    SIZE(isendblkstart,2),SIZE(isendblksize,2),SIZE(ioffset,2))

      !-------------------------------------------------------------------------
      !  Determine intersecting mesh blocks
      !-------------------------------------------------------------------------
      dosend = .TRUE.
      DO k=1,pdim
          fromistartk=this%istart(k,isub)+offset(k)
          fromnnodesk=this%nnodes(k,isub)-chop(k)
          toistartk  =to_mesh%istart(k,jsub)-ghostsize(k)
          tonnodesk  =to_mesh%nnodes(k,jsub)+2*ghostsize(k)
          iblockstart(k) = MAX(fromistartk,toistartk)
          iblockstopk = MIN((fromistartk+fromnnodesk),(toistartk+tonnodesk))
          nblocksize(k)  = iblockstopk-iblockstart(k)
          ! do not send if there is nothing to be sent
          IF (nblocksize(k) .LT. 1) dosend = .FALSE.
      ENDDO
      IF (dosend) THEN
          nsendlist = nsendlist + 1
          IF (nsendlist .GT. isize) THEN
              !-----------------------------------------------------------------
              !  Grow memory for the sendlists
              !-----------------------------------------------------------------
              isize  = 2*isize
              iopt   = ppm_param_alloc_grow_preserve
              ldu(1) = isize
              CALL ppm_alloc(isendfromsub,ldu,iopt,info)
                    or_fail_alloc("isendfromsub")
              CALL ppm_alloc(isendtosub,ldu,iopt,info)
                    or_fail_alloc("isendtosub")
              ldu(1) = pdim
              ldu(2) = isize
              CALL ppm_alloc(isendblkstart,ldu,iopt,info)
                    or_fail_alloc("isendblkstart")
              CALL ppm_alloc(isendblksize,ldu,iopt,info)
                    or_fail_alloc("isendblksize")
              CALL ppm_alloc(ioffset,ldu,iopt,info)
                    or_fail_alloc("ioffset")
          ENDIF
          !---------------------------------------------------------------------
          !  send block (isendblkstart...isendblkstart+nblocksize-1) 
          !  from sub isub to sub jsub
          !---------------------------------------------------------------------
          ! global sub index of local sub where the blocks come from
          isendfromsub(nsendlist)        = isub
          ! global sub index of other sub where the blocks go to
          isendtosub(nsendlist)          = jsub
          ! start of mesh block in global mesh
          isendblkstart(1:pdim,nsendlist)= iblockstart(1:pdim) - offset(1:pdim)
          ! size of mesh block
          isendblksize(1:pdim,nsendlist) = nblocksize(1:pdim)
          ! mesh block offset for periodic images
          ioffset(1:pdim,nsendlist)      = offset(1:pdim)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
      end_subroutine()

      RETURN
      CONTAINS
      SUBROUTINE check
          IF (SIZE(ghostsize,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'ghostsize must be at least of length ppm_dim',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (SIZE(offset,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'offset must be at least of length ppm_dim',__LINE__,info)
              GOTO 8888
          ENDIF
          IF ((isub.LE.0).OR.(isub.GT.ppm_topo(this%topoID)%t%nsubs)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'isub out of range',__LINE__,info)
              GOTO 8888
          ENDIF
          IF ((jsub.LE.0).OR.(jsub.GT.ppm_topo(to_mesh%topoID)%t%nsubs)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'jsub out of range',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE equi_mesh_block_intersect

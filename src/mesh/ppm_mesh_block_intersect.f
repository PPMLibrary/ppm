      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_mesh_block_intersect
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
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

      SUBROUTINE ppm_mesh_block_intersect(from_topoid,to_topoid,from_meshid,   &
     &    to_meshid,isub,jsub,offset,ghostsize,nsendlist,isendfromsub,         &
     &    isendtosub,isendblkstart,isendblksize,ioffset,info,lsymm)
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
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_typedef
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_check_id
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: from_topoid
      !!! Source topology identifier
      INTEGER                 , INTENT(IN   ) :: to_topoid
      !!! Target topology identifier
      INTEGER                 , INTENT(IN   ) :: to_meshid
      !!! Target mesh identifier
      INTEGER                 , INTENT(IN   ) :: from_meshid
      !!! Source mesh identifier
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
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: valid
      TYPE(ppm_t_equi_mesh), POINTER   :: from_mesh => NULL()
      TYPE(ppm_t_equi_mesh), POINTER   :: to_mesh   => NULL()
      TYPE(ppm_t_topo),      POINTER   :: from_topo => NULL()
      TYPE(ppm_t_topo),      POINTER   :: to_topo   => NULL()
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_mesh_block_intersect',t0,info)
      pdim = ppm_dim

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      from_topo => ppm_topo(from_topoid)%t
      to_topo => ppm_topo(to_topoid)%t
      from_mesh => from_topo%mesh(from_meshid)
      to_mesh => to_topo%mesh(to_meshid)

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
          fromistartk=from_mesh%istart(k,isub)+offset(k)
          fromnnodesk=from_mesh%nnodes(k,isub)-chop(k)
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
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_mesh_block_intersect',   &
     &                'local send source sub list ISENDFROMSUB',__LINE__,info)
                  GOTO 9999
              ENDIF
              CALL ppm_alloc(isendtosub,ldu,iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_mesh_block_intersect',   &
     &               'local send destination sub list ISENDTOSUB',__LINE__,info)
                  GOTO 9999
              ENDIF
              ldu(1) = pdim
              ldu(2) = isize
              CALL ppm_alloc(isendblkstart,ldu,iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_mesh_block_intersect',   &
     &               'local send block start list ISENDBLKSTART',__LINE__,info)
                  GOTO 9999
              ENDIF
              CALL ppm_alloc(isendblksize,ldu,iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_mesh_block_intersect',   &
     &               'local send block size list ISENDBLKSIZE',__LINE__,info)
                  GOTO 9999
              ENDIF
              CALL ppm_alloc(ioffset,ldu,iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_mesh_block_intersect',   &
     &               'local send block offset list IOFFSET',__LINE__,info)
                  GOTO 9999
              ENDIF
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
 9999 CONTINUE
      CALL substop('ppm_mesh_block_intersect',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        CALL ppm_check_topoid(from_topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'from_topoid out of range',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(to_topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'to_topoid out of range',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_meshid(from_topoid,from_meshid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'from_meshid out of range',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_meshid(to_topoid,to_meshid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'to_meshid out of range',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (SIZE(ghostsize,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'ghostsize must be at least of length ppm_dim',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (SIZE(offset,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'offset must be at least of length ppm_dim',__LINE__,info)
              GOTO 8888
          ENDIF
          IF ((isub.LE.0).OR.(isub.GT.ppm_topo(from_topoid)%t%nsubs)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'isub out of range',__LINE__,info)
              GOTO 8888
          ENDIF
          IF ((jsub.LE.0).OR.(jsub.GT.ppm_topo(to_topoid)%t%nsubs)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'jsub out of range',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_mesh_block_intersect

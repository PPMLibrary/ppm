      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_field_init
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

      SUBROUTINE ppm_map_field_init(topoid,target_topoid,meshid, &
     &                              target_meshid,info)
      !!! This routine must be called before the ppm_map_field_globalstored
      !!! routine.
      !!!
      !!! [WARNING]
      !!! This routine has not been tested, reviewed, or checked. Comments
      !!! and documentation are wrong. This routine might kill your cat or
      !!! worse.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_typedef
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_check_id
      USE ppm_module_mesh_block_intersect
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! Topology ID of source
      !!!
      !!! CAUTION: used to be target
      INTEGER                 , INTENT(IN   ) :: target_topoid
      !!! Topology ID of target
      INTEGER                 , INTENT(IN   ) :: meshid
      !!! Mesh ID of source
      INTEGER                 , INTENT(IN   ) :: target_meshid
      !!! Mesh ID of target
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)            :: ldu
      INTEGER, DIMENSION(ppm_dim)      :: op
      INTEGER, DIMENSION(ppm_dim,27)   :: ond
      INTEGER, DIMENSION(ppm_dim)      :: iblockstart,nblocksize
      INTEGER, DIMENSION(ppm_dim)      :: ghostsize
      INTEGER                          :: i,j,k,idom,sendrank,recvrank
      INTEGER                          :: iopt,iset,ibuffer,pdim,nnd
      INTEGER                          :: nsendlist,nsend
      INTEGER                          :: nrecvlist
      CHARACTER(ppm_char)              :: mesg
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: lsouth,lnorth,least,lwest,ltop,lbottom
      LOGICAL, DIMENSION(3)            :: lsymm
      LOGICAL                          :: valid
      TYPE(ppm_t_topo),      POINTER   :: topo
      TYPE(ppm_t_topo),      POINTER   :: target_topo
      TYPE(ppm_t_equi_mesh), POINTER   :: mesh
      TYPE(ppm_t_equi_mesh), POINTER   :: target_mesh
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_field_init',t0,info)
      pdim = ppm_dim

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t
      target_topo => ppm_topo(target_topoid)%t
      mesh => topo%mesh(meshid)
      target_mesh => target_topo%mesh(meshid)


      IF (ppm_buffer_set .GT. 0) THEN
        info = ppm_error_warning
        CALL ppm_error(ppm_err_map_incomp,'ppm_map_field_init',  &
     &      'Buffer was not empty. Possible loss of data!',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls (used for checks!)
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_global

      !-------------------------------------------------------------------------
      !  Check if origin and target meshes are compatible (i.e. have the
      !  same number of grid points in the whole comput. domain)
      !-------------------------------------------------------------------------
      DO i=1,pdim
          IF (mesh%Nm(i) .NE. target_mesh%Nm(i)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_bad_mesh,'ppm_map_field_init',  &
     &            'source and destination meshes are incompatible',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDDO


      !-------------------------------------------------------------------------
      !  Pre-calculate shift offsets for periodic ghost images of subs
      !-------------------------------------------------------------------------
      op(1) = mesh%Nm(1) - 1
      op(2) = mesh%Nm(2) - 1
      IF (pdim .GT. 2) THEN
          op(3) = mesh%Nm(3) - 1
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the local temporary sendlists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = topo%nsublist
      CALL ppm_alloc(isendfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'local send source sub list ISENDFROMSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(isendtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'local send destination sub list ISENDTOSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = pdim
      ldu(2) = topo%nsublist
      CALL ppm_alloc(isendblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'local send block start list ISENDBLKSTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(isendblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'local send block size list ISENDBLKSIZE',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ioffset,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'local send block offset list IOFFSET',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find intersecting mesh domains to be sent
      !-------------------------------------------------------------------------
      nsendlist = 0
      ghostsize = 0
      DO i=1,topo%nsublist
          idom = topo%isublist(i)
          ond(1:pdim,1:27) = 0
          !---------------------------------------------------------------------
          !  We now use symmetry. Only send up to ndata-1
          !  unless we are at non-periodic, non-internal boundary.
          !  Because we do this, we need to consider periodic images.
          !  In periodic systems, the box needs to be shifted and overlaps
          !  recomputed in order to get the periodic image neighbors.
          !  Check if  face of this sub coincides with a periodic domain
          !  boundary.
          !---------------------------------------------------------------------
          lsymm(1:pdim) = .TRUE.
          IF ((topo%bcdef(2) .NE. ppm_param_bcdef_periodic) .AND. &
     &        (topo%subs_bc(2,idom) .NE. 0)) lsymm(1) = .FALSE.
          IF ((topo%bcdef(4) .NE. ppm_param_bcdef_periodic) .AND. &
     &        (topo%subs_bc(4,idom) .NE. 0)) lsymm(2) = .FALSE.
          IF ((topo%bcdef(6) .NE. ppm_param_bcdef_periodic) .AND. &
     &        (topo%subs_bc(6,idom) .NE. 0)) lsymm(3) = .FALSE.

          lwest  = .FALSE.
          least  = .FALSE.
          IF ((topo%bcdef(1) .EQ. ppm_param_bcdef_periodic) .AND. &
     &        (topo%subs_bc(1,idom) .NE. 0)) lwest = .TRUE.
          lsouth = .FALSE.
          lnorth = .FALSE.
          IF ((topo%bcdef(3) .EQ. ppm_param_bcdef_periodic) .AND. &
     &        (topo%subs_bc(3,idom) .NE. 0)) lsouth = .TRUE.
          IF (pdim .GT. 2) THEN
              lbottom= .FALSE.
              ltop   = .FALSE.
              IF ((topo%bcdef(5) .EQ. ppm_param_bcdef_periodic) .AND. &
     &            (topo%subs_bc(5,idom) .NE. 0)) lbottom = .TRUE.
          ENDIF

          !---------------------------------------------------------------------
          !  Determine number of shifts and actual shift indices needed
          !---------------------------------------------------------------------
          nnd = 1
          IF (lwest) THEN
              nnd = nnd + 1
              ond(1,nnd) = op(1)
          ENDIF
          IF (least) THEN
              nnd = nnd + 1
              ond(1,nnd) = -op(1)
          ENDIF
          IF (lsouth) THEN
              nnd = nnd + 1
              ond(2,nnd) = op(2)
              IF (lwest) THEN
                  nnd = nnd + 1
                  ond(1,nnd) = op(1)
                  ond(2,nnd) = op(2)
              ENDIF
              IF (least) THEN
                  nnd = nnd + 1
                  ond(1,nnd) = -op(1)
                  ond(2,nnd) = op(2)
              ENDIF
          ENDIF
          IF (lnorth) THEN
              nnd = nnd + 1
              ond(2,nnd) = -op(2)
              IF (lwest) THEN
                  nnd = nnd + 1
                  ond(1,nnd) = op(1)
                  ond(2,nnd) = -op(2)
              ENDIF
              IF (least) THEN
                  nnd = nnd + 1
                  ond(1,nnd) = -op(1)
                  ond(2,nnd) = -op(2)
              ENDIF
          ENDIF
          IF (pdim .GT. 2) THEN
              IF (lbottom) THEN
                  nnd = nnd + 1
                  ond(3,nnd) = op(3)
                  IF (lwest) THEN
                      nnd = nnd + 1
                      ond(1,nnd) = op(1)
                      ond(3,nnd) = op(3)
                  ENDIF
                  IF (least) THEN
                      nnd = nnd + 1
                      ond(1,nnd) = -op(1)
                      ond(3,nnd) = op(3)
                  ENDIF
                  IF (lsouth) THEN
                      nnd = nnd + 1
                      ond(2,nnd) = op(2)
                      ond(3,nnd) = op(3)
                      IF (lwest) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = op(1)
                          ond(2,nnd) = op(2)
                          ond(3,nnd) = op(3)
                      ENDIF
                      IF (least) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = -op(1)
                          ond(2,nnd) = op(2)
                          ond(3,nnd) = op(3)
                      ENDIF
                  ENDIF
                  IF (lnorth) THEN
                      nnd = nnd + 1
                      ond(2,nnd) = -op(2)
                      ond(3,nnd) = op(3)
                      IF (lwest) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = op(1)
                          ond(2,nnd) = -op(2)
                          ond(3,nnd) = op(3)
                      ENDIF
                      IF (least) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = -op(1)
                          ond(2,nnd) = -op(2)
                          ond(3,nnd) = op(3)
                      ENDIF
                  ENDIF
              ENDIF
              IF (ltop) THEN
                  nnd = nnd + 1
                  ond(3,nnd) = -op(3)
                  IF (lwest) THEN
                      nnd = nnd + 1
                      ond(1,nnd) = op(1)
                      ond(3,nnd) = -op(3)
                  ENDIF
                  IF (least) THEN
                      nnd = nnd + 1
                      ond(1,nnd) = -op(1)
                      ond(3,nnd) = -op(3)
                  ENDIF
                  IF (lsouth) THEN
                      nnd = nnd + 1
                      ond(2,nnd) = op(2)
                      ond(3,nnd) = -op(3)
                      IF (lwest) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = op(1)
                          ond(2,nnd) = op(2)
                          ond(3,nnd) = -op(3)
                      ENDIF
                      IF (least) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = -op(1)
                          ond(2,nnd) = op(2)
                          ond(3,nnd) = -op(3)
                      ENDIF
                  ENDIF
                  IF (lnorth) THEN
                      nnd = nnd + 1
                      ond(2,nnd) = -op(2)
                      ond(3,nnd) = -op(3)
                      IF (lwest) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = op(1)
                          ond(2,nnd) = -op(2)
                          ond(3,nnd) = -op(3)
                      ENDIF
                      IF (least) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = -op(1)
                          ond(2,nnd) = -op(2)
                          ond(3,nnd) = -op(3)
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF

          !---------------------------------------------------------------------
          !  Do the shifting and compute mesh blocks
          !---------------------------------------------------------------------
          DO k=1,nnd
              DO j=1,target_topo%nsubs
                  CALL ppm_mesh_block_intersect(topoid,target_topoid,meshid,&
     &                target_meshid,idom,j,ond(1:pdim,k),ghostsize,nsendlist, &
     &                isendfromsub,isendtosub,isendblkstart,isendblksize, &
     &                ioffset,info,lsymm)
                  IF (info .NE. 0) GOTO 9999
              ENDDO
          ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for the local temporary receive lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = target_topo%nsublist
      CALL ppm_alloc(irecvfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'local recv source sub list IRECVFROMSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(irecvtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'local recv destination sub list IRECVTOSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = pdim
      ldu(2) = target_topo%nsublist
      CALL ppm_alloc(irecvblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'local recv block start list IRECVBLKSTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(irecvblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'local recv block size list IRECVBLKSIZE',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find intersecting mesh domains to be received
      !-------------------------------------------------------------------------
      nrecvlist = 0
      ! loop over fromtopo first and THEN over totopo in order to get the
      ! same ordering of mesh blocks as in the sendlist. This is crutial
      ! for the push and the pop to work properly !!!
      DO j=1,topo%nsubs
          ond(1:pdim,1:27) = 0
          !---------------------------------------------------------------------
          !  We now use symmetry. Only send up to ndata-1
          !  unless we are at non-periodic, non-internal boundary.
          !  Because we do this, we need to consider periodic images.
          !  In periodic systems, the box needs to be shifted and overlaps
          !  recomputed in order to get the periodic image neighbors.
          !  Check if  face of this sub coincides with a periodic domain
          !  boundary.
          !---------------------------------------------------------------------
          lsymm(1:pdim) = .TRUE.
          IF ((topo%bcdef(2) .NE. ppm_param_bcdef_periodic) .AND. &
     &        (topo%subs_bc(2,j) .NE. 0)) lsymm(1) = .FALSE.
          IF ((topo%bcdef(4) .NE. ppm_param_bcdef_periodic) .AND. &
     &        (topo%subs_bc(4,j) .NE. 0)) lsymm(2) = .FALSE.
          IF ((topo%bcdef(6) .NE. ppm_param_bcdef_periodic) .AND. &
     &        (topo%subs_bc(6,j) .NE. 0)) lsymm(3) = .FALSE.

          lwest  = .FALSE.
          least  = .FALSE.
          IF ((topo%bcdef(1) .EQ. ppm_param_bcdef_periodic) .AND. &
     &        (topo%subs_bc(1,j) .NE. 0)) lwest = .TRUE.
          lsouth = .FALSE.
          lnorth = .FALSE.
          IF ((topo%bcdef(3) .EQ. ppm_param_bcdef_periodic) .AND. &
     &        (topo%subs_bc(3,j) .NE. 0)) lsouth = .TRUE.
          IF (pdim .GT. 2) THEN
              lbottom= .FALSE.
              ltop   = .FALSE.
              IF ((topo%bcdef(5) .EQ. ppm_param_bcdef_periodic) .AND. &
     &            (topo%subs_bc(5,j) .NE. 0)) lbottom = .TRUE.
          ENDIF

          !---------------------------------------------------------------------
          !  Determine number of shifts and actual shift indices needed
          !---------------------------------------------------------------------
          nnd = 1
          IF (lwest) THEN
              nnd = nnd + 1
              ond(1,nnd) = op(1)
          ENDIF
          IF (least) THEN
              nnd = nnd + 1
              ond(1,nnd) = -op(1)
          ENDIF
          IF (lsouth) THEN
              nnd = nnd + 1
              ond(2,nnd) = op(2)
              IF (lwest) THEN
                  nnd = nnd + 1
                  ond(1,nnd) = op(1)
                  ond(2,nnd) = op(2)
              ENDIF
              IF (least) THEN
                  nnd = nnd + 1
                  ond(1,nnd) = -op(1)
                  ond(2,nnd) = op(2)
              ENDIF
          ENDIF
          IF (lnorth) THEN
              nnd = nnd + 1
              ond(2,nnd) = -op(2)
              IF (lwest) THEN
                  nnd = nnd + 1
                  ond(1,nnd) = op(1)
                  ond(2,nnd) = -op(2)
              ENDIF
              IF (least) THEN
                  nnd = nnd + 1
                  ond(1,nnd) = -op(1)
                  ond(2,nnd) = -op(2)
              ENDIF
          ENDIF
          IF (pdim .GT. 2) THEN
              IF (lbottom) THEN
                  nnd = nnd + 1
                  ond(3,nnd) = op(3)
                  IF (lwest) THEN
                      nnd = nnd + 1
                      ond(1,nnd) = op(1)
                      ond(3,nnd) = op(3)
                  ENDIF
                  IF (least) THEN
                      nnd = nnd + 1
                      ond(1,nnd) = -op(1)
                      ond(3,nnd) = op(3)
                  ENDIF
                  IF (lsouth) THEN
                      nnd = nnd + 1
                      ond(2,nnd) = op(2)
                      ond(3,nnd) = op(3)
                      IF (lwest) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = op(1)
                          ond(2,nnd) = op(2)
                          ond(3,nnd) = op(3)
                      ENDIF
                      IF (least) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = -op(1)
                          ond(2,nnd) = op(2)
                          ond(3,nnd) = op(3)
                      ENDIF
                  ENDIF
                  IF (lnorth) THEN
                      nnd = nnd + 1
                      ond(2,nnd) = -op(2)
                      ond(3,nnd) = op(3)
                      IF (lwest) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = op(1)
                          ond(2,nnd) = -op(2)
                          ond(3,nnd) = op(3)
                      ENDIF
                      IF (least) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = -op(1)
                          ond(2,nnd) = -op(2)
                          ond(3,nnd) = op(3)
                      ENDIF
                  ENDIF
              ENDIF
              IF (ltop) THEN
                  nnd = nnd + 1
                  ond(3,nnd) = -op(3)
                  IF (lwest) THEN
                      nnd = nnd + 1
                      ond(1,nnd) = op(1)
                      ond(3,nnd) = -op(3)
                  ENDIF
                  IF (least) THEN
                      nnd = nnd + 1
                      ond(1,nnd) = -op(1)
                      ond(3,nnd) = -op(3)
                  ENDIF
                  IF (lsouth) THEN
                      nnd = nnd + 1
                      ond(2,nnd) = op(2)
                      ond(3,nnd) = -op(3)
                      IF (lwest) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = op(1)
                          ond(2,nnd) = op(2)
                          ond(3,nnd) = -op(3)
                      ENDIF
                      IF (least) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = -op(1)
                          ond(2,nnd) = op(2)
                          ond(3,nnd) = -op(3)
                      ENDIF
                  ENDIF
                  IF (lnorth) THEN
                      nnd = nnd + 1
                      ond(2,nnd) = -op(2)
                      ond(3,nnd) = -op(3)
                      IF (lwest) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = op(1)
                          ond(2,nnd) = -op(2)
                          ond(3,nnd) = -op(3)
                      ENDIF
                      IF (least) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = -op(1)
                          ond(2,nnd) = -op(2)
                          ond(3,nnd) = -op(3)
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF

          !---------------------------------------------------------------------
          !  Do the shifting and compute mesh blocks
          !---------------------------------------------------------------------
          DO k=1,nnd
              DO i=1,target_topo%nsublist
                  idom = target_topo%isublist(i)
                  CALL ppm_mesh_block_intersect(topoid,target_topoid,meshid,&
     &                target_meshid,j,idom,ond(1:pdim,k),ghostsize,nrecvlist, &
     &                irecvfromsub,irecvtosub,irecvblkstart,irecvblksize, &
     &                ioffset,info,lsymm)
                 IF (info .NE. 0) GOTO 9999
              ENDDO
          ENDDO
      ENDDO


      !-------------------------------------------------------------------------
      !  Allocate memory for the global mesh sendlists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldu(1) = nsendlist
      CALL ppm_alloc(mesh%mapping%isendfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'source send sub list PPM_MESH_ISENDFROMSUB',__LINE__,info)
          GOTO 9999
      ENDIF

      ldu(1) = pdim
      ldu(2) = nsendlist
      CALL ppm_alloc(mesh%mapping%isendblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'send block start list PPM_MESH_ISENDBLKSTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(mesh%mapping%isendblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'send block size list PPM_MESH_ISENDBLKSIZE',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the global mesh receive lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldu(1) = nrecvlist
      CALL ppm_alloc(mesh%mapping%irecvtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'destination recv sub list PPM_MESH_IRECVTOSUB',__LINE__,info)
          GOTO 9999
      ENDIF

      ldu(1) = pdim
      ldu(2) = nrecvlist
      CALL ppm_alloc(mesh%mapping%irecvblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'recv block start list PPM_MESH_IRECVBLKSTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(mesh%mapping%irecvblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'recv block size list PPM_MESH_IRECVBLKSIZE',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the global send/recv lists
      !-------------------------------------------------------------------------
      ldu(1) = ppm_nproc
      CALL ppm_alloc(mesh%mapping%isendlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'global send rank list PPM_ISENDLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(mesh%mapping%irecvlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'global recv rank list PPM_IRECVLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = ppm_nproc + 1
      CALL ppm_alloc(mesh%mapping%psendbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'global send buffer pointer PPM_PSENDBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(mesh%mapping%precvbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_init',     &
     &        'global recv buffer pointer PPM_PRECVBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Reset the number of buffer entries
      !-------------------------------------------------------------------------
      ppm_buffer_set = 0

      !-------------------------------------------------------------------------
      !  loop over all processors; First is the processor itself
      !-------------------------------------------------------------------------
      sendrank           = ppm_rank - 1
      recvrank           = ppm_rank + 1
      mesh%mapping%psendbuffer(1) = 1
      mesh%mapping%precvbuffer(1) = 1
      mesh%mapping%nsendlist      = 0
      mesh%mapping%nrecvlist      = 0
      iset               = 0
      ibuffer            = 0

      DO i=1,ppm_nproc
         !----------------------------------------------------------------------
         !  compute the next processor
         !----------------------------------------------------------------------
         sendrank = sendrank + 1
         IF (sendrank.GE.ppm_nproc) sendrank = sendrank - ppm_nproc
         recvrank = recvrank - 1
         IF (recvrank.LT.        0) recvrank = recvrank + ppm_nproc

         !----------------------------------------------------------------------
         !  Store the processor to which we will send to
         !----------------------------------------------------------------------
         mesh%mapping%nsendlist = mesh%mapping%nsendlist + 1
         mesh%mapping%isendlist(mesh%mapping%nsendlist) = sendrank

         !----------------------------------------------------------------------
         !  Store the processor to which we will recv from
         !----------------------------------------------------------------------
         mesh%mapping%nrecvlist = mesh%mapping%nrecvlist + 1
         mesh%mapping%irecvlist(mesh%mapping%nrecvlist) = recvrank

         !----------------------------------------------------------------------
         !  Find all mesh blocks that are sent and store them.
         !  To get all mesh blocks that are sent in round
         !  i=1,ppm_nsendlist:
         !     dest_rank = ppm_isendlist(i)
         !     ilo = ppm_psendbuffer(i)
         !     ihi = ppm_psendbuffer(i+1)-1
         !     DO j=ilo,ihi
         !       SEND(ppm_mesh_isendblk(:,j) of sub ppm_mesh_isendfromsub(j) to
         !         processor dest_rank to sub ppm_mesh_isendtosub(j))
         !     ENDDO
         !----------------------------------------------------------------------
         ibuffer = mesh%mapping%nsendlist + 1
         mesh%mapping%psendbuffer(ibuffer) = &
     &                     mesh%mapping%psendbuffer(mesh%mapping%nsendlist)
         DO j=1,nsendlist
             IF (target_topo%sub2proc(isendtosub(j)) .EQ. sendrank) THEN
                 mesh%mapping%psendbuffer(ibuffer) = &
     &                               mesh%mapping%psendbuffer(ibuffer) + 1
                 iset = mesh%mapping%psendbuffer(ibuffer) - 1
                 mesh%mapping%isendfromsub(iset) = isendfromsub(j)
                 mesh%mapping%isendblkstart(1:pdim,iset)=isendblkstart(1:pdim,j)
                 mesh%mapping%isendblksize(1:pdim,iset)  =isendblksize(1:pdim,j)
                 IF (ppm_debug .GT. 1) THEN
                     IF (ppm_dim .EQ. 2) THEN
                         WRITE(mesg,'(2(A,2I4),A,I3)') ' sending ',  &
     &                   isendblkstart(1:2,j),' of size ',isendblksize(1:2,j),&
     &                   ' to ',sendrank
                     ELSEIF (ppm_dim .EQ. 3) THEN
                         WRITE(mesg,'(2(A,3I4),A,I3)') ' sending ',  &
     &                   isendblkstart(1:3,j),' of size ',isendblksize(1:3,j),&
     &                   ' to ',sendrank
                     ENDIF
                     CALL ppm_write(ppm_rank,'ppm_map_field_init',mesg,info)
                 ENDIF
             ENDIF
         ENDDO

         !----------------------------------------------------------------------
         !  Find all mesh blocks that are received and store them.
         !  To get all mesh blocks that are received in round
         !  i=1,ppm_nrecvlist:
         !     source_rank = ppm_irecvlist(i)
         !     ilo = ppm_precvbuffer(i)
         !     ihi = ppm_precvbuffer(i+1)-1
         !     DO j=ilo,ihi
         !       RECV(ppm_mesh_irecvblk(:,j) of sub ppm_mesh_irecvfromsub(j)
         !         from processor source_rank to sub ppm_mesh_irecvtosub(j))
         !     ENDDO
         !----------------------------------------------------------------------
         ibuffer = mesh%mapping%nrecvlist + 1
         mesh%mapping%precvbuffer(ibuffer) = &
     &                     mesh%mapping%precvbuffer(mesh%mapping%nrecvlist)
         DO j=1,nrecvlist
             IF (topo%sub2proc(irecvfromsub(j)) .EQ. recvrank) THEN
                 mesh%mapping%precvbuffer(ibuffer) = &
     &                                 mesh%mapping%precvbuffer(ibuffer) + 1
                 iset = mesh%mapping%precvbuffer(ibuffer) - 1
                 mesh%mapping%irecvtosub(iset) = irecvtosub(j)
                 mesh%mapping%irecvblkstart(1:pdim,iset) = &
     &                               irecvblkstart(1:pdim,j) + ioffset(1:pdim,j)
                 mesh%mapping%irecvblksize(1:pdim,iset) = irecvblksize(1:pdim,j)
                 IF (ppm_debug .GT. 1) THEN
                     IF (ppm_dim .EQ. 2) THEN
                         WRITE(mesg,'(2(A,2I4),A,I3)') ' recving ',  &
     &                   irecvblkstart(1:2,j),' of size ',irecvblksize(1:2,j),&
     &                   ' from ',recvrank
                     ELSEIF (ppm_dim .EQ. 3) THEN
                         WRITE(mesg,'(2(A,3I4),A,I3)') ' recving ',  &
     &                   irecvblkstart(1:3,j),' of size ',irecvblksize(1:3,j),&
     &                   ' from ',recvrank
                     ENDIF
                     CALL ppm_write(ppm_rank,'ppm_map_field_init',mesg,info)
                 ENDIF
             ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Deallocate memory of the local lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      CALL ppm_alloc(isendfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_init',     &
     &        'local send source sub list ISENDFROMSUB',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isendtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_init',     &
     &        'local send destination sub list ISENDTOSUB',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isendblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_init',     &
     &        'local send block start list ISENDBLKSTART',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isendblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_init',     &
     &        'local send block size list ISENDBLKSIZE',__LINE__,info)
      ENDIF
      CALL ppm_alloc(irecvfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_init',     &
     &        'local recv source sub list IRECVFROMSUB',__LINE__,info)
      ENDIF
      CALL ppm_alloc(irecvtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_init',     &
     &        'local recv destination sub list IRECVTOSUB',__LINE__,info)
      ENDIF
      CALL ppm_alloc(irecvblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_init',     &
     &        'local recv block start list IRECVBLKSTART',__LINE__,info)
      ENDIF
      CALL ppm_alloc(irecvblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_init',     &
     &        'local recv block size list IRECVBLKSIZE',__LINE__,info)
      ENDIF
      CALL ppm_alloc(ioffset,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_init',     &
     &        'local recv block offset list IOFFSET',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_field_init',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          CALL ppm_check_topoid(topoid,ltop,info)
          IF (.NOT. ltop) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_ghost_init',  &
     &            'topoid out of range',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_meshid(topoid,meshid,ltop,info)
          IF (.NOT. ltop) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_ghost_init',  &
     &            'meshid out of range',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(target_topoid,ltop,info)
          IF (.NOT. ltop) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_ghost_init',  &
     &            'topoid out of range',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_meshid(target_topoid,target_meshid,ltop,info)
          IF (.NOT. ltop) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_ghost_init',  &
     &            'meshid out of range',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (SIZE(ghostsize,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_ghost_init',  &
     &            'ghostsize must be at least of length ppm_dim',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_map_field_init

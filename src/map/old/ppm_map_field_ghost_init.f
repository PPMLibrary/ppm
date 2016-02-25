      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_map_field_ghost_init
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

      SUBROUTINE ppm_map_field_ghost_init(topoid,meshid,ghostsize,info)
      !!! This routine sets up the lists needed for sending and receiving
      !!! ghost mesh points on a given topology/mesh combination.
      !!! These lists are stored for later use by ppm_map_field_ghost_get and
      !!! `ppm_map_field_ghost_put`.
      !!!
      !!! [IMPORTANT]
      !!! This routine should never be called by the user, the library will
      !!! call this routine before the first `ppm_map_field_ghost_get` or
      !!! `ppm_map_field_ghost_put`.
      !!!
      !!! [NOTE]
      !!! The lists are now saved in the mesh structure for which the ghost
      !!! mappings are initialized. There is no need to reinitialize after
      !!! doing global mappings
      !!!
      !!! [WARNING]
      !!! The new version of this routine might cause some memory overhead
      !!!
      !!! [NOTE]
      !!! The first part of the send/recv lists contains the on-processor data.
      !!! The names send and recv are correct for the ghost get map. For
      !!! the put case, they need to be interpreted the other way round.
      !!! The shift cases for periodic images are currently hard coded.
      !!! Probably there is a more elegant way of doing this.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_mpi
      USE ppm_module_mesh_block_intersect
      USE ppm_module_util_commopt
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! Topology ID
      INTEGER                 , INTENT(IN   ) :: meshid
      !!! Mesh ID
      INTEGER, DIMENSION(:)   , INTENT(IN   ) :: ghostsize
      !!! Size of the ghost layer in numbers of grid points in all space
      !!! dimensions (1...ppm_dim).
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(4)            :: ldu
      INTEGER, DIMENSION(ppm_dim)      :: op
      INTEGER, DIMENSION(ppm_dim,26)   :: ond
      INTEGER                          :: i,j,sendrank,recvrank,idom,jdom,k
      INTEGER                          :: iopt,iset,ibuffer,pdim,isize,nnd
      INTEGER                          :: nsendlist,nsend,tag1,lb,ub,nrecv
      CHARACTER(ppm_char)              :: mesg
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: lsouth,lnorth,least,lwest,ltop,lbottom
#ifdef __MPI
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: commstat
#endif
      TYPE(ppm_t_topo), POINTER        :: topo
      TYPE(ppm_t_equi_mesh), POINTER   :: mesh
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_field_ghost_init',t0,info)
      pdim = ppm_dim

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t

      SELECT TYPE (t => ppm_mesh%vec(meshid)%t)
      TYPE IS (ppm_t_equi_mesh)
         mesh => t
      END SELECT

      !-------------------------------------------------------------------------
      !  check if the optimal communication protocol is known for this
      !  topology
      !-------------------------------------------------------------------------
      IF (.NOT.topo%isoptimized) THEN
        !-----------------------------------------------------------------------
        !  if not: determine it before calling map_field_ghost_init
        !-----------------------------------------------------------------------
        CALL ppm_util_commopt(topoid,info)
        IF (info.NE.0) GOTO 9999
        IF (ppm_debug .GT. 1) THEN
            DO i=1,topo%nneighproc
                WRITE(mesg,'(A,I4)') 'have neighbor: ',topo%ineighproc(i)
                CALL ppm_write(ppm_rank,'ppm_map_field_ghost_init',mesg,info)
            END DO
            DO i=1,topo%ncommseq
                WRITE(mesg,'(A,I4)') 'communicate: ',topo%icommseq(i)
                CALL ppm_write(ppm_rank,'ppm_map_field_ghost_init',mesg,info)
            END DO
        ENDIF
      ENDIF


      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls (used to check push/pop)
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_init

      !-------------------------------------------------------------------------
      !  Allocate memory for the local temporary sendlists
      !-------------------------------------------------------------------------
      isize  = topo%nsublist
      iopt   = ppm_param_alloc_fit
      ldu(1) = isize
      CALL ppm_alloc(isendfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'local send source sub list ISENDFROMSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(isendtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'local send destination sub list ISENDTOSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = pdim
      ldu(2) = isize
      CALL ppm_alloc(isendblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'local send block start list ISENDBLKSTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(isendblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'local send block size list ISENDBLKSIZE',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ioffset,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'local mesh block offset list IOFFSET',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Pre-calculate shift offsets for periodic ghost images of subs
      !-------------------------------------------------------------------------

      op(1) = mesh%Nm(1) - 1
      op(2) = mesh%Nm(2) - 1
      IF (pdim .GT. 2) THEN
          op(3) = mesh%Nm(3) - 1
      ENDIF

      !-------------------------------------------------------------------------
      !  Find intersecting mesh domains to be sent for the ghost get
      !-------------------------------------------------------------------------
      nsendlist = 0
      DO i=1,topo%nsublist
          idom = topo%isublist(i)
          ond(1:pdim,1:26) = 0
          !---------------------------------------------------------------------
          !  First, we check the neighbors of the subdomain
          !---------------------------------------------------------------------
          DO j=1,topo%nneighsubs(i)
              jdom = topo%ineighsubs(j,i)
              ! source and destination meshes and topologies are identical
              CALL ppm_mesh_block_intersect(topoid,topoid,meshid,meshid,       &
     &           idom,jdom,ond(1:pdim,1),ghostsize,nsendlist,isendfromsub,     &
     &           isendtosub,isendblkstart,isendblksize,ioffset,info)
              IF (info .NE. 0) GOTO 9999
          ENDDO

          !---------------------------------------------------------------------
          !  In periodic systems, the box needs to be shifted and overlaps
          !  recomputed in order to get the periodic image neighbors.
          !  Check if any face of this sub coincides with a periodic domain
          !  boundary.
          !---------------------------------------------------------------------
          lwest  = .FALSE.
          IF ((topo%bcdef(1) .EQ. ppm_param_bcdef_periodic) .AND. &
     &        (topo%subs_bc(1,idom) .NE. 0)) lwest = .TRUE.
          least  = .FALSE.
          IF ((topo%bcdef(2) .EQ. ppm_param_bcdef_periodic) .AND. &
     &        (topo%subs_bc(2,idom) .NE. 0)) least = .TRUE.
          lsouth = .FALSE.
          IF ((topo%bcdef(3) .EQ. ppm_param_bcdef_periodic) .AND. &
     &        (topo%subs_bc(3,idom) .NE. 0)) lsouth = .TRUE.
          lnorth = .FALSE.
          IF ((topo%bcdef(4) .EQ. ppm_param_bcdef_periodic) .AND. &
     &        (topo%subs_bc(4,idom) .NE. 0)) lnorth = .TRUE.
          IF (pdim .GT. 2) THEN
              lbottom= .FALSE.
              IF ((topo%bcdef(5) .EQ. ppm_param_bcdef_periodic) .AND. &
     &            (topo%subs_bc(5,idom) .NE. 0)) lbottom = .TRUE.
              ltop   = .FALSE.
              IF ((topo%bcdef(6) .EQ. ppm_param_bcdef_periodic) .AND. &
     &            (topo%subs_bc(6,idom) .NE. 0)) ltop = .TRUE.
          ENDIF

          !---------------------------------------------------------------------
          !  Determine number of shifts and actual shift indices needed
          !---------------------------------------------------------------------
          nnd = 0
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
          !  Do the shifting and compute periodic ghost mesh blocks
          !---------------------------------------------------------------------
          DO k=1,nnd
              ! first with the original (non-shifted) image of itself
              jdom = idom
              CALL ppm_mesh_block_intersect(topoid,topoid,meshid,meshid,       &
     &           idom,jdom,ond(1:pdim,k),ghostsize,nsendlist,isendfromsub,     &
     &           isendtosub,isendblkstart,isendblksize,ioffset,info)
              IF (info .NE. 0) GOTO 9999
              ! Then with all the neighbors
              DO j=1,topo%nneighsubs(i)
                  jdom = topo%ineighsubs(j,i)
                  CALL ppm_mesh_block_intersect(topoid,topoid,meshid,meshid,   &
     &                idom,jdom,ond(1:pdim,k),ghostsize,nsendlist,isendfromsub,&
     &                isendtosub,isendblkstart,isendblksize,ioffset,info)
                  IF (info .NE. 0) GOTO 9999
              ENDDO
          ENDDO
      ENDDO    ! i=1,ppm_nsublist

      !-------------------------------------------------------------------------
      !  Grow memory for the mesh ghost maps
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldu(1) = nsendlist
      CALL ppm_alloc(mesh%ghost_fromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost source subs MESH%GHOST_FROMSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      !CALL ppm_alloc(ppm_mesh_ghost_recvfromsub,ldu,iopt,info)
      !IF (info .NE. 0) THEN
      !    info = ppm_error_fatal
      !    CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
   !  &        'ghost recv source subs PPM_MESH_GHOST_RECVFROMSUB',__LINE__,info)
   !       GOTO 9999
   !   ENDIF
      CALL ppm_alloc(mesh%ghost_tosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost target subs MESH%GHOST_TOSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(mesh%ghost_recvtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost recv target subs MESH%GHOST_RECVTOSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = topo%ncommseq + 1
      CALL ppm_alloc(mesh%ghost_blk,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost block list PMESH%GHOST_BLK',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(mesh%ghost_recvblk,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost recv block list MESH%GHOST_RECVBLK',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = pdim
      ldu(2) = nsendlist
      CALL ppm_alloc(mesh%ghost_blkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost block start list MESH%GHOST_BLKSTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(mesh%ghost_recvblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost recv block start list MESH%GHOST_RECVBLKSTART',  &
     &        __LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(mesh%ghost_blksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost block size list MESH%GHOST_BLKSIZE',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(mesh%ghost_recvblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost recv block size list MESH%GHOST_RECVBLKSIZE', &
     &        __LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate local memory for sorted offset list
      !-------------------------------------------------------------------------
      ldu(1) = pdim
      ldu(2) = nsendlist
      CALL ppm_alloc(mesh_ghost_offset,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost block offset list MESH_GHOST_OFFSET',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  loop over the neighboring processors according to the optimized
      !  communication sequence.
      !-------------------------------------------------------------------------
      mesh%ghost_nsend = nsendlist
      mesh%ghost_nrecv = 0
      mesh%ghost_blk(1) = 1
      mesh%ghost_recvblk(1) = 1
      iset               = 0
      ibuffer            = 0
      isize              = SIZE(mesh%ghost_fromsub,1)

      DO i=1,topo%ncommseq
         !----------------------------------------------------------------------
         !  get the next processor to send/recv with
         !----------------------------------------------------------------------
         sendrank = topo%icommseq(i)
         recvrank = sendrank

         !----------------------------------------------------------------------
         !  Initialize mesh block pointers
         !----------------------------------------------------------------------
         ibuffer = i + 1
         mesh%ghost_blk(ibuffer) = mesh%ghost_blk(i)
         mesh%ghost_recvblk(ibuffer) = mesh%ghost_recvblk(i)

         !----------------------------------------------------------------------
         !  Only assign data if there is any communication in this round
         !----------------------------------------------------------------------
         IF (sendrank .GE. 0) THEN

             !------------------------------------------------------------------
             !  Find all mesh blocks that are sent and store them.
             !------------------------------------------------------------------
             DO j=1,nsendlist
                 IF (topo%sub2proc(isendtosub(j)) .EQ. sendrank) THEN
                     mesh%ghost_blk(ibuffer) = mesh%ghost_blk(ibuffer)+1
                     iset = mesh%ghost_blk(ibuffer) - 1
                     ! store this for the topology as it can be reused
                     mesh%ghost_fromsub(iset)=isendfromsub(j)
                     mesh%ghost_tosub(iset)  =isendtosub(j)
                     mesh%ghost_blkstart(1:pdim,iset)=isendblkstart(1:pdim,j)
                     mesh%ghost_blksize(1:pdim,iset) =isendblksize(1:pdim,j)
                     ! also re-order the offsets as we need them for
                     ! computing the receive lists further down !!
                     mesh_ghost_offset(1:pdim,iset) = ioffset(1:pdim,j)
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
                         CALL ppm_write(ppm_rank,'ppm_map_field_ghost_init',  &
     &                       mesg,info)
                     ENDIF
                 ENDIF
             ENDDO

             !------------------------------------------------------------------
             !  Build receive lists for the on-processor part of the data
             !------------------------------------------------------------------
             IF (sendrank .EQ. ppm_rank) THEN
                 ub = mesh%ghost_blk(2)
                 mesh%ghost_nrecv = ub - 1
                 mesh%ghost_recvblk(2) = ub
                 DO j=1,ub-1
                     mesh%ghost_recvtosub(j) = mesh%ghost_tosub(j)
  !                   ppm_mesh_ghost_recvfromsub(j,meshid,topoid) =  &
  !   &                   ppm_mesh_ghost_fromsub(j,meshid,topoid)
                     mesh%ghost_recvblkstart(1:pdim,j) = &
     &                   mesh%ghost_blkstart(1:pdim,j) + &
     &                   mesh_ghost_offset(1:pdim,j)
                     mesh%ghost_recvblksize(1:pdim,j) = &
     &                   mesh%ghost_blksize(1:pdim,j)
                 ENDDO

#ifdef __MPI
             ELSE
                 !--------------------------------------------------------------
                 !  Communicate the block indices
                 !--------------------------------------------------------------
                 lb = mesh%ghost_blk(i)
                 ub = mesh%ghost_blk(ibuffer)
                 ! How many blocks am I sending to that guy
                 nsend = ub - lb
                 tag1 = 100
                 CALL MPI_SendRecv(nsend,1,MPI_INTEGER,sendrank,tag1,nrecv,1,  &
     &               MPI_INTEGER,recvrank,tag1,ppm_comm,commstat,info)
                 ! How many blocks will I receive from the guy?
                 mesh%ghost_nrecv = mesh%ghost_nrecv + nrecv
                 mesh%ghost_recvblk(ibuffer) = mesh%ghost_recvblk(i) + nrecv

                 !--------------------------------------------------------------
                 !  Check if receive lists are large enough
                 !--------------------------------------------------------------
                 IF (mesh%ghost_nrecv .GT. isize) THEN
                     !----------------------------------------------------------
                     !  Grow receive lists if needed
                     !----------------------------------------------------------
                     isize = mesh%ghost_nrecv
                     iopt = ppm_param_alloc_grow_preserve
                     ldu(1) = isize
   !                  CALL ppm_alloc(ppm_mesh_ghost_recvfromsub,ldu,iopt,info)
   !                  IF (info .NE. 0) THEN
   !                      info = ppm_error_fatal
   !                      CALL ppm_error(ppm_err_alloc,   &
   !  &                     'ppm_map_field_ghost_init',     &
   !  &                     'ghost recv source subs PPM_MESH_GHOST_RECVFROMSUB',&
   !  &                     __LINE__,info)
   !                      GOTO 9999
   !                  ENDIF
                     CALL ppm_alloc(mesh%ghost_recvtosub,ldu,iopt,info)
                     IF (info .NE. 0) THEN
                         info = ppm_error_fatal
                         CALL ppm_error(ppm_err_alloc,   &
     &                     'ppm_map_field_ghost_init',     &
     &                     'ghost recv target subs PPM_MESH_GHOST_RECVTOSUB',&
     &                     __LINE__,info)
                         GOTO 9999
                     ENDIF
                     ldu(1) = pdim
                     ldu(2) = isize
                     CALL ppm_alloc(mesh%ghost_recvblkstart,ldu,iopt,info)
                     IF (info .NE. 0) THEN
                         info = ppm_error_fatal
                         CALL ppm_error(ppm_err_alloc,   &
     &                    'ppm_map_field_ghost_init',     &
     &                    'ghost recv block start PPM_MESH_GHOST_RECVBLKSTART',&
     &                     __LINE__,info)
                         GOTO 9999
                     ENDIF
                     CALL ppm_alloc(mesh%ghost_recvblksize,ldu,iopt,info)
                     IF (info .NE. 0) THEN
                         info = ppm_error_fatal
                         CALL ppm_error(ppm_err_alloc,   &
     &                    'ppm_map_field_ghost_init',     &
     &                    'ghost recv block size PPM_MESH_GHOST_RECVBLKSIZE',&
     &                     __LINE__,info)
                         GOTO 9999
                     ENDIF
                 ENDIF

                 !--------------------------------------------------------------
                 !  Allocate memory for block data send and recv buffers
                 !--------------------------------------------------------------
                 iopt   = ppm_param_alloc_grow
                 ldu(1) = nsend*(2*pdim+1)
                 CALL ppm_alloc(sendbuf,ldu,iopt,info)
                 IF (info .NE. 0) THEN
                     info = ppm_error_fatal
                     CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init', &
     &                   'local MPI send buffer SENDBUF',__LINE__,info)
                     GOTO 9999
                 ENDIF
                 ldu(1) = nrecv*(2*pdim+1)
                 CALL ppm_alloc(recvbuf,ldu,iopt,info)
                 IF (info .NE. 0) THEN
                     info = ppm_error_fatal
                     CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init', &
     &                   'local MPI recv buffer RECVBUF',__LINE__,info)
                     GOTO 9999
                 ENDIF

                 !--------------------------------------------------------------
                 !  Pack and send all the ghost mesh block data
                 !--------------------------------------------------------------
                 ! Pack all the send data
                 iset = 0
                 DO j=lb,ub-1
                     iset = iset + 1
                     sendbuf(iset) = mesh%ghost_tosub(j)
             !        iset = iset + 1
             !        sendbuf(iset) = ppm_mesh_ghost_fromsub(j,meshid,topoid)
                     sendbuf(iset+1:iset+pdim) =    &
     &                   mesh%ghost_blkstart(1:pdim,j) +  &
     &                   mesh_ghost_offset(1:pdim,j)
                     iset = iset + pdim
                     sendbuf(iset+1:iset+pdim) =    &
     &                   mesh%ghost_blksize(1:pdim,j)
                     iset = iset + pdim
                 ENDDO
                 ! Send it to the destination processor and get my stuff
                 tag1 = 200
                 CALL MPI_SendRecv(sendbuf,iset,MPI_INTEGER,sendrank,tag1, &
     &                             recvbuf,nrecv*(2*pdim+1),MPI_INTEGER,   &
     &                             recvrank,tag1,ppm_comm,commstat,info)
                 ! Unpack the received data
                 lb = mesh%ghost_recvblk(i)
                 ub = mesh%ghost_recvblk(ibuffer)
                 iset = 0
                 DO j=lb,ub-1
                     iset = iset + 1
                     mesh%ghost_recvtosub(j) = recvbuf(iset)
                 !    iset = iset + 1
                 !    ppm_mesh_ghost_recvfromsub(j,meshid,topoid) = recvbuf(iset)
                     mesh%ghost_recvblkstart(1:pdim,j) =   &
     &                   recvbuf(iset+1:iset+pdim)
                     iset = iset + pdim
                     mesh%ghost_recvblksize(1:pdim,j) =    &
     &                   recvbuf(iset+1:iset+pdim)
                     iset = iset + pdim

                 ENDDO
#endif
             ENDIF ! sendrank.EQ.ppm_rank
          ENDIF    ! sendrank .GE. 0
      ENDDO

      !-------------------------------------------------------------------------
      !  Deallocate memory of the local lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      CALL ppm_alloc(recvbuf,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'local MPI recv buffer RECVBUF',__LINE__,info)
      ENDIF
      CALL ppm_alloc(sendbuf,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'local MPI send buffer SENDBUF',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isendfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'local send source sub list ISENDFROMSUB',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isendtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'local send destination sub list ISENDTOSUB',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isendblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'local send block start list ISENDBLKSTART',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isendblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'local send block size list ISENDBLKSIZE',__LINE__,info)
      ENDIF
      CALL ppm_alloc(ioffset,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'local send block offset list IOFFSET',__LINE__,info)
      ENDIF
      CALL ppm_alloc(mesh_ghost_offset,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'ghost block offset list MESH_GHOST_OFFSET',__LINE__,info)
      ENDIF


      mesh%ghost_initialized = .TRUE.

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_field_ghost_init',t0,info)
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
          IF (SIZE(ghostsize,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_ghost_init',  &
     &            'ghostsize must be at least of length ppm_dim',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_map_field_ghost_init

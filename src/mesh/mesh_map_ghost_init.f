      SUBROUTINE equi_mesh_map_ghost_init(this,info,ghostsize)
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
      USE ppm_module_mpi
      USE ppm_module_util_commopt
      USE ppm_module_topo_typedef,    ONLY : ppm_topo,ppm_t_topo
      USE ppm_module_interfaces, ONLY : ppm_t_mesh_mapping_
      USE ppm_module_mapping_typedef, ONLY : ppm_t_mesh_mapping
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_t_equi_mesh)                         :: this

      INTEGER,                         INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: ghostsize
      !!! Size of the ghost layer in numbers of grid points in all space
      !!! dimensions (1...ppm_dim).
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), POINTER :: topo

      CLASS(ppm_t_mesh_mapping_), POINTER :: map

      INTEGER, DIMENSION(2)          :: ldu
      INTEGER, DIMENSION(ppm_dim)    :: op
      INTEGER, DIMENSION(ppm_dim)    :: ghostsize_
      INTEGER, DIMENSION(ppm_dim,26) :: ond
      INTEGER                        :: i,j,sendrank,recvrank,isub,jsub,k
      INTEGER                        :: iopt,iset,ibuffer,isize,nnd
      INTEGER                        :: nsendlist,nsend,tag,lb,ub,nrecv
#ifdef __MPI
      INTEGER                        :: requestr
      INTEGER, DIMENSION(2)          :: requests
#endif

      LOGICAL :: lsouth,lnorth,least,lwest,ltop,lbottom

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      start_subroutine("mesh_map_ghost_init")

      topo => ppm_topo(this%topoid)%t

      !-------------------------------------------------------------------------
      !  check if the optimal communication protocol is known for this
      !  topology
      !-------------------------------------------------------------------------
      IF (.NOT.topo%isoptimized) THEN
         !-----------------------------------------------------------------------
         !  if not: determine it before calling map_field_ghost_init
         !-----------------------------------------------------------------------
         CALL ppm_util_commopt(this%topoid,info)
         IF (info.NE.0) GOTO 9999
         IF (ppm_debug.GT.1) THEN
            DO i=1,topo%nneighproc
               stdout_f('(A,I4)',"have neighbor: ",'topo%ineighproc(i)')
            ENDDO
            DO i=1,topo%ncommseq
               stdout_f('(A,I4)',"communicate: ",'topo%icommseq(i)')
            ENDDO
         ENDIF
      ENDIF

      IF (PRESENT(ghostsize)) THEN
         ghostsize_=MIN(ghostsize,this%ghostsize)
      ELSE
         ghostsize_=this%ghostsize
      ENDIF

      ALLOCATE(ppm_t_mesh_mapping::map,STAT=info)
      or_fail_alloc("Failed to allocate map")

      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_dim
      CALL ppm_alloc(map%ghostsize,ldu,iopt,info)
      or_fail_alloc("map%ghostsize")

      map%ghostsize=ghostsize_

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
      or_fail_alloc("isendfromsub")

      CALL ppm_alloc(isendtosub,ldu,iopt,info)
      or_fail_alloc("isendtosub")

      ldu(1) = ppm_dim
      ldu(2) = isize
      CALL ppm_alloc(isendpatchid,ldu,iopt,info)
      or_fail_alloc("isendpatchid")

      CALL ppm_alloc(isendblkstart,ldu,iopt,info)
      or_fail_alloc("isendblkstart")

      CALL ppm_alloc(isendblksize,ldu,iopt,info)
      or_fail_alloc("isendblksize")

      CALL ppm_alloc(ioffset,ldu,iopt,info)
      or_fail_alloc("ioffset")

      !-------------------------------------------------------------------------
      !  Pre-calculate shift offsets for periodic ghost images of subs
      !-------------------------------------------------------------------------
      op=this%Nm-1

      !-------------------------------------------------------------------------
      !  Find intersecting mesh domains to be sent for the ghost get
      !-------------------------------------------------------------------------
      nsendlist = 0
      DO i=1,topo%nsublist
         isub = topo%isublist(i)
         ond  = 0
         !---------------------------------------------------------------------
         !  First, we check the neighbors of the subdomain
         !---------------------------------------------------------------------
         DO j=1,topo%nneighsubs(i)
            jsub = topo%ineighsubs(j,i)
            ! source and destination meshes and topologies are identical
            CALL this%block_intersect(this,isub,jsub,ond(1:ppm_dim,1), &
            &    nsendlist,isendfromsub,isendtosub,isendpatchid,    &
            &    isendblkstart,isendblksize,ioffset,info,ghostsize=ghostsize_)
            or_fail("block_intersect failed")
         ENDDO

         !---------------------------------------------------------------------
         !  In periodic systems, the box needs to be shifted and overlaps
         !  recomputed in order to get the periodic image neighbors.
         !  Check if any face of this sub coincides with a periodic domain
         !  boundary.
         !---------------------------------------------------------------------
         lwest  = .FALSE.
         IF ((topo%bcdef(1) .EQ. ppm_param_bcdef_periodic) .AND. &
         &   (topo%subs_bc(1,isub) .NE. 0)) lwest = .TRUE.
         least  = .FALSE.
         IF ((topo%bcdef(2) .EQ. ppm_param_bcdef_periodic) .AND. &
         &   (topo%subs_bc(2,isub) .NE. 0)) least = .TRUE.
         lsouth = .FALSE.
         IF ((topo%bcdef(3) .EQ. ppm_param_bcdef_periodic) .AND. &
         &   (topo%subs_bc(3,isub) .NE. 0)) lsouth = .TRUE.
         lnorth = .FALSE.
         IF ((topo%bcdef(4) .EQ. ppm_param_bcdef_periodic) .AND. &
         &   (topo%subs_bc(4,isub) .NE. 0)) lnorth = .TRUE.
         IF (ppm_dim .GT. 2) THEN
            lbottom= .FALSE.
            IF ((topo%bcdef(5) .EQ. ppm_param_bcdef_periodic) .AND. &
            &   (topo%subs_bc(5,isub) .NE. 0)) lbottom = .TRUE.
            ltop   = .FALSE.
            IF ((topo%bcdef(6) .EQ. ppm_param_bcdef_periodic) .AND. &
            &   (topo%subs_bc(6,isub) .NE. 0)) ltop = .TRUE.
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
         IF (ppm_dim .GT. 2) THEN
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
            CALL this%block_intersect(this,isub,isub,ond(1:ppm_dim,k), &
            &    nsendlist,isendfromsub,isendtosub,isendpatchid,    &
            &    isendblkstart,isendblksize,ioffset,info,ghostsize=ghostsize_)
            or_fail("block_intersect failed")
            ! Then with all the neighbors
            DO j=1,topo%nneighsubs(i)
               jsub = topo%ineighsubs(j,i)
               CALL this%block_intersect(this,isub,jsub,ond(1:ppm_dim,k), &
               &    nsendlist,isendfromsub,isendtosub,isendpatchid,    &
               &    isendblkstart,isendblksize,ioffset,info,ghostsize=ghostsize_)
               or_fail("block_intersect failed")
            ENDDO
         ENDDO
      ENDDO    ! i=1,ppm_nsublist

      !-------------------------------------------------------------------------
      !  Grow memory for the mesh ghost maps
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldu(1) = nsendlist
      CALL ppm_alloc(map%ghost_fromsub,ldu,iopt,info)
      or_fail_alloc("map%ghost_fromsub")

      CALL ppm_alloc(map%ghost_tosub,ldu,iopt,info)
      or_fail_alloc("map%ghost_tosub")

      CALL ppm_alloc(map%ghost_recvtosub,ldu,iopt,info)
      or_fail_alloc("map%ghost_recvtosub")

      ldu(1) = topo%ncommseq + 1
      CALL ppm_alloc(map%ghost_blk,ldu,iopt,info)
      or_fail_alloc("map%ghost_blk")

      CALL ppm_alloc(map%ghost_recvblk,ldu,iopt,info)
      or_fail_alloc("map%ghost_recvblk")

      ldu(1) = ppm_dim
      ldu(2) = nsendlist
      CALL ppm_alloc(map%ghost_patchid,ldu,iopt,info)
      or_fail_alloc("map%ghost_patchid")

      CALL ppm_alloc(map%ghost_recvpatchid,ldu,iopt,info)
      or_fail_alloc("map%recvghost_patchid")

      CALL ppm_alloc(map%ghost_blkstart,ldu,iopt,info)
      or_fail_alloc("map%ghost_blkstart")

      CALL ppm_alloc(map%ghost_recvblkstart,ldu,iopt,info)
      or_fail_alloc("map%ghost_recvblkstart")

      CALL ppm_alloc(map%ghost_blksize,ldu,iopt,info)
      or_fail_alloc("map%ghost_blksize")

      CALL ppm_alloc(map%ghost_recvblksize,ldu,iopt,info)
      or_fail_alloc("map%ghost_recvblksize")

      !-------------------------------------------------------------------------
      !  Allocate local memory for sorted offset list
      !-------------------------------------------------------------------------
      ldu(1) = ppm_dim
      ldu(2) = nsendlist
      CALL ppm_alloc(mesh_ghost_offset,ldu,iopt,info)
      or_fail_alloc("mesh_ghost_offset,ldu")

      !-------------------------------------------------------------------------
      !  loop over the neighboring processors according to the optimized
      !  communication sequence.
      !-------------------------------------------------------------------------
      map%ghost_nsend      = nsendlist
      map%ghost_nrecv      = 0
      map%ghost_blk(1)     = 1
      map%ghost_recvblk(1) = 1
      iset                 = 0
      ibuffer              = 0
      isize                = SIZE(map%ghost_fromsub,1)

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

         map%ghost_blk(ibuffer)     = map%ghost_blk(i)
         map%ghost_recvblk(ibuffer) = map%ghost_recvblk(i)

         !----------------------------------------------------------------------
         !  Only assign data if there is any communication in this round
         !----------------------------------------------------------------------
         IF (sendrank .GE. 0) THEN

            !------------------------------------------------------------------
            !  Find all mesh blocks that are sent and store them.
            !------------------------------------------------------------------
            DO j=1,nsendlist
               IF (topo%sub2proc(isendtosub(j)) .EQ. sendrank) THEN
                  map%ghost_blk(ibuffer) = map%ghost_blk(ibuffer)+1
                  iset = map%ghost_blk(ibuffer) - 1
                  ! store this for the topology as it can be reused
                  map%ghost_fromsub(iset)         = isendfromsub(j)
                  map%ghost_tosub(iset)           = isendtosub(j)
                  map%ghost_patchid(1:ppm_dim,iset)  = isendpatchid(1:ppm_dim,j)
                  map%ghost_blkstart(1:ppm_dim,iset) = isendblkstart(1:ppm_dim,j)
                  map%ghost_blksize(1:ppm_dim,iset)  = isendblksize(1:ppm_dim,j)
                  ! also re-order the offsets as we need them for
                  ! computing the receive lists further down !!
                  mesh_ghost_offset(1:ppm_dim,iset)  = ioffset(1:ppm_dim,j)
                  IF (ppm_debug .GT. 1) THEN
                     IF (ppm_dim .EQ. 2) THEN
                        stdout_f('(2(A,2I4),3(A,I0))'," sending ",'isendblkstart(1:2,j)', &
                        & " of size ",'isendblksize(1:2,j)'," on sub ",'isendfromsub(j)', &
                        & " to sub ",'isendtosub(j)'," on proc ",sendrank)
                     ELSEIF (ppm_dim .EQ. 3) THEN
                        stdout_f('(2(A,3I4),3(A,I0))'," sending ",'isendblkstart(1:3,j)', &
                        & " of size ",'isendblksize(1:3,j)'," on sub ",'isendfromsub(j)', &
                        & " to sub ",'isendtosub(j)'," on proc ",sendrank)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

            !------------------------------------------------------------------
            !  Build receive lists for the on-processor part of the data
            !------------------------------------------------------------------
            IF (sendrank .EQ. ppm_rank) THEN
               ub = map%ghost_blk(2)
               map%ghost_nrecv = ub - 1
               map%ghost_recvblk(2) = ub
               DO j=1,ub-1
                  map%ghost_recvtosub(j)           = map%ghost_tosub(j)
                  map%ghost_recvpatchid(1:ppm_dim,j)  = map%ghost_patchid(1:ppm_dim,j)
                  map%ghost_recvblkstart(1:ppm_dim,j) = map%ghost_blkstart(1:ppm_dim,j) &
                  &                                + mesh_ghost_offset(1:ppm_dim,j)
                  map%ghost_recvblksize(1:ppm_dim,j)  = map%ghost_blksize(1:ppm_dim,j)
               ENDDO
#ifdef __MPI
            ELSE
               IF (ppm_rank.GT.recvrank) THEN
                  tag=ppm_rank*(ppm_rank-1)/2+recvrank+100
               ELSE
                  tag=recvrank*(recvrank-1)/2+ppm_rank+100
               ENDIF

               CALL MPI_Irecv(nrecv,1,MPI_INTEGER,recvrank,tag,ppm_comm,requestr,info)
               or_fail_MPI("MPI_Irecv")

               !--------------------------------------------------------------
               !  Communicate the block indices
               !--------------------------------------------------------------
               lb = map%ghost_blk(i)
               ub = map%ghost_blk(ibuffer)
               ! How many blocks am I sending to that guy
               nsend = ub - lb

               IF (ppm_rank.GT.sendrank) THEN
                  tag=ppm_rank*(ppm_rank-1)/2+sendrank+100
               ELSE
                  tag=sendrank*(sendrank-1)/2+ppm_rank+100
               ENDIF

               CALL MPI_Isend(nsend,1,MPI_INTEGER,sendrank,tag,ppm_comm,requests(1),info)
               or_fail_MPI("MPI_Isend")

               !--------------------------------------------------------------
               !  Allocate memory for block data send buffers
               !--------------------------------------------------------------
               iopt   = ppm_param_alloc_grow
               ldu(1) = nsend*(3*ppm_dim+1)
               CALL ppm_alloc(sendbuf,ldu,iopt,info)
               or_fail_alloc("sendbuf")

               !--------------------------------------------------------------
               !  Pack and send all the ghost mesh block data
               !--------------------------------------------------------------
               ! Pack all the send data
               iset = 0
               DO j=lb,ub-1
                  iset = iset + 1
                  sendbuf(iset) = map%ghost_tosub(j)
                  sendbuf(iset+1:iset+ppm_dim) = map%ghost_patchid(1:ppm_dim,j)

                  iset = iset + ppm_dim
                  sendbuf(iset+1:iset+ppm_dim) = map%ghost_blkstart(1:ppm_dim,j) &
                  &                         + mesh_ghost_offset(1:ppm_dim,j)

                  iset = iset + ppm_dim
                  sendbuf(iset+1:iset+ppm_dim) = map%ghost_blksize(1:ppm_dim,j)

                  iset = iset + ppm_dim
               ENDDO

               ! Send it to the destination processor and get my stuff
               IF (ppm_rank.GT.recvrank) THEN
                  tag=ppm_rank*(ppm_rank-1)/2+recvrank+200
               ELSE
                  tag=recvrank*(recvrank-1)/2+ppm_rank+200
               ENDIF

               !Wait to have nrecv
               CALL MPI_Wait(requestr,MPI_STATUS_IGNORE,info)
               or_fail_MPI("MPI_Wait")

               !--------------------------------------------------------------
               !  Allocate memory for block data send and recv buffers
               !--------------------------------------------------------------
               iopt   = ppm_param_alloc_grow
               ldu(1) = nrecv*(3*ppm_dim+1)
               CALL ppm_alloc(recvbuf,ldu,iopt,info)
               or_fail_alloc("recvbuf")

               ! Send it to the destination processor and get my stuff
               CALL MPI_Irecv(recvbuf,ldu(1),MPI_INTEGER,recvrank,tag,ppm_comm,requestr,info)
               or_fail_MPI("MPI_Irecv")

               IF (ppm_rank.GT.sendrank) THEN
                  tag=ppm_rank*(ppm_rank-1)/2+sendrank+200
               ELSE
                  tag=sendrank*(sendrank-1)/2+ppm_rank+200
               ENDIF

               CALL MPI_Isend(sendbuf,iset,MPI_INTEGER,sendrank,tag,ppm_comm,requests(2),info)
               or_fail_MPI("MPI_Isend")

               ! How many blocks will I receive from the guy?
               map%ghost_nrecv            = map%ghost_nrecv      + nrecv
               map%ghost_recvblk(ibuffer) = map%ghost_recvblk(i) + nrecv
               !--------------------------------------------------------------
               !  Check if receive lists are large enough
               !--------------------------------------------------------------
               IF (map%ghost_nrecv .GT. isize) THEN
                  !----------------------------------------------------------
                  !  Grow receive lists if needed
                  !----------------------------------------------------------
                  isize = map%ghost_nrecv
                  iopt = ppm_param_alloc_grow_preserve
                  ldu(1) = isize
                  CALL ppm_alloc(map%ghost_recvtosub,ldu,iopt,info)
                  or_fail_alloc("map%ghost_recvtosub")
                  ldu(1) = ppm_dim
                  ldu(2) = isize
                  CALL ppm_alloc(map%ghost_recvblkstart,ldu,iopt,info)
                  or_fail_alloc("map%ghost_recvblkstart")

                  CALL ppm_alloc(map%ghost_recvblksize,ldu,iopt,info)
                  or_fail_alloc("map%ghost_recvblksize")

                  CALL ppm_alloc(map%ghost_recvpatchid,ldu,iopt,info)
                  or_fail_alloc("map%recvghost_patchid")
               ENDIF

               ! Unpack the received data
               lb = map%ghost_recvblk(i)
               ub = map%ghost_recvblk(ibuffer)

               CALL MPI_Wait(requestr,MPI_STATUS_IGNORE,info)
               or_fail_MPI("MPI_Wait")

               iset = 0
               DO j=lb,ub-1
                  iset = iset + 1
                  map%ghost_recvtosub(j)           = recvbuf(iset)
                  map%ghost_recvpatchid(1:ppm_dim,j)  = recvbuf(iset+1:iset+ppm_dim)

                  iset = iset + ppm_dim
                  map%ghost_recvblkstart(1:ppm_dim,j) = recvbuf(iset+1:iset+ppm_dim)

                  iset = iset + ppm_dim
                  map%ghost_recvblksize(1:ppm_dim,j)  = recvbuf(iset+1:iset+ppm_dim)

                  iset = iset + ppm_dim
               ENDDO

               CALL MPI_Waitall(2,requests,MPI_STATUSES_IGNORE,info)
               or_fail_MPI("MPI_Waitall")
#endif
             ENDIF ! sendrank.EQ.ppm_rank
         ENDIF ! sendrank .GE. 0
      ENDDO !i=1,topo%ncommseq

      !-------------------------------------------------------------------------
      !  Deallocate memory of the local lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      CALL ppm_alloc(recvbuf,ldu,iopt,info)
      or_fail_dealloc("recvbuf")

      CALL ppm_alloc(sendbuf,ldu,iopt,info)
      or_fail_dealloc("sendbuf")

      CALL ppm_alloc(isendfromsub,ldu,iopt,info)
      or_fail_dealloc("isendfromsub")

      CALL ppm_alloc(isendtosub,ldu,iopt,info)
      or_fail_dealloc("isendtosub")

      CALL ppm_alloc(isendpatchid,ldu,iopt,info)
      or_fail_dealloc("isendpatchid")

      CALL ppm_alloc(isendblkstart,ldu,iopt,info)
      or_fail_dealloc("isendblkstart")

      CALL ppm_alloc(isendblksize,ldu,iopt,info)
      or_fail_dealloc("isendblksize")

      CALL ppm_alloc(ioffset,ldu,iopt,info)
      or_fail_dealloc("ioffset")

      CALL ppm_alloc(mesh_ghost_offset,ldu,iopt,info)
      or_fail_dealloc("mesh_ghost_offset")

      CALL this%maps%push(map,info)
      or_fail("Failed to push the map!")

      end_subroutine()

      END SUBROUTINE equi_mesh_map_ghost_init

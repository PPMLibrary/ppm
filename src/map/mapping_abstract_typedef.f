#if __MYTYPE == __MAPPINGTYPE
      !----------------------------------------------------------------------
      !  Generic mapping data type
      !----------------------------------------------------------------------
      TYPE,ABSTRACT :: ppm_t_mapping_
          INTEGER                           :: source_topoid = -1
          !!! topology ID on which the data was stored before computing the mapping
          INTEGER                           :: target_topoid = -1
          !!! topology ID on which the data should be stored after the mapping
          !------------------------------------------------------------------
          !  buffers for communication
          !------------------------------------------------------------------
          INTEGER,  DIMENSION(:),   POINTER :: nsend => NULL()
          INTEGER,  DIMENSION(:),   POINTER :: nrecv => NULL()
          INTEGER,  DIMENSION(:),   POINTER :: psend => NULL()
          INTEGER,  DIMENSION(:),   POINTER :: precv => NULL()
          INTEGER,  DIMENSION(:,:), POINTER :: pp    => NULL()
          INTEGER,  DIMENSION(:,:), POINTER :: qq    => NULL()

          INTEGER                           :: old_nsendlist = 0
          INTEGER                           :: old_buffer_set = 0

          !------------------------------------------------------------------
          !
          !------------------------------------------------------------------
          INTEGER,  DIMENSION(:),   POINTER :: ppm_psendbuffer => NULL()
          !!! pointer to particles within the send buffer
          !!!
          !!! In terms of particle *not* the actual position in the buffer
          INTEGER,  DIMENSION(:),   POINTER :: ppm_precvbuffer => NULL()
          !!! pointer to particles within the recv buffer
          !!!
          !!! In terms of particle *not* the actual position in the buffer
          INTEGER                           :: ppm_nsendbuffer
          !!! the size of the send buffer (not the array but the number of
          !!! used elements)
          !!!
          !!! In terms of entries in the buffer *not* the number of particles
          INTEGER                           :: ppm_nrecvbuffer
          !!! the size of the recv buffer (not the array but the number of
          !!! used elements)
          !!!
          !!! In terms of entries in the buffer *not* the number of particles
          INTEGER                           :: ppm_sendbufsize
          !!! the actual size of the send buffer array
          INTEGER                           :: ppm_recvbufsize
          !!! the actual size of the receive buffer array
          INTEGER                           :: ppm_buffer_set
          !!! the total number of particle fields packed in
          !!! the send buffer, ie. xp, vp is two sets
          INTEGER,  DIMENSION(:),   POINTER :: ppm_buffer_type  => NULL()
          !!! types of the data on the sendbuffer
          INTEGER,  DIMENSION(:),   POINTER :: ppm_buffer_dim   => NULL()
          !!! dimensions of the original on-processor particle arrays.
          !----------------------------------------------------------------------
          !  mapping
          !----------------------------------------------------------------------
          INTEGER                           :: ppm_map_type
          INTEGER                           :: ppm_nsendlist
          INTEGER                           :: ppm_nrecvlist
          INTEGER,  DIMENSION(:),   POINTER :: ppm_isendlist => NULL()
          INTEGER,  DIMENSION(:),   POINTER :: ppm_irecvlist => NULL()
      END TYPE ppm_t_mapping_
#endif

      TYPE,ABSTRACT,EXTENDS(ppm_t_mapping_) :: DTYPE(ppm_t_part_mapping)_
          !------------------------------------------------------------------
          !  buffers for communication
          !------------------------------------------------------------------
          REAL(MK), DIMENSION(:),   POINTER :: send  => NULL()
          REAL(MK), DIMENSION(:),   POINTER :: recv  => NULL()

          INTEGER                           :: oldNpart
          INTEGER                           :: newNpart

          INTEGER,  DIMENSION(:,:), POINTER :: ppm_ghosthack    => NULL()
          !!! invert map of ghost for symmetry
          INTEGER,  DIMENSION(:),   POINTER :: ppm_buffer2part  => NULL()
          !!! Used for the original on-processor particle IDs in the order in which
          !!! they are in the sendbuffer. Used to push additional particle
          !!! data on the buffer in the correct order.

          !------------------------------------------------------------------
          !
          !------------------------------------------------------------------
          REAL(MK), DIMENSION(:),   POINTER :: ppm_recvbuffer  => NULL()
          !!! receive buffer for particles
          REAL(MK), DIMENSION(:),   POINTER :: ppm_sendbuffer  => NULL()
          !!! send buffer for particles
          REAL(MK), DIMENSION(:),   POINTER :: ppm_ghost_offset => NULL()
          !!! ghost offset
          !!!
          !!! ghost particles may have a spatial offset compared to their real
          !!! particle - we need to store this offset in terms of the buffer id
          !!! in order to be able to push/send/pop ghost coordinates in the
          !!! ghost_get and ghost_put mapping; this is needed when using Verlet
          !!! lists: here nothing is remapped and for symmetry we need both to
          !!! put the particle forces etc. and to get the updated particle
          !!! position; for the asymmetric case we do not put the forces back but
          !!! we need to get the updated positions of the ghost - NOT getting
          !!! ghosts - because a new mapping is NOT needed but simply reusing
          !!! the old mapping push/send/popping the updated coordinates
          !!! The offsets are stored as the pdata in the psendbuffer, i.e.,
          !!! `ppm_ghost_offset(ibuffer+0) = xp_offset(1,)`                     +
          !!! `ppm_ghost_offset(ibuffer+1) = xp_offset(2,)`                     +
          !!! `ppm_ghost_offset(ibuffer+2) = xp_offset(3,)`
          REAL(MK), DIMENSION(:),   POINTER :: ppm_ghost_offset_fac => NULL()
          !!! ghost offset factor
          !!!
          !!! This buffer is needed for the special case of symmetric and
          !!! antisymmetric boundary conditions, because the offset is not a
          !!! simple constant to be added to the particle coordinate
      CONTAINS
          PROCEDURE(DTYPE(map_part_create)_),  DEFERRED :: create
          PROCEDURE(DTYPE(map_part_destroy)_), DEFERRED :: destroy

      END TYPE DTYPE(ppm_t_part_mapping)_
minclude ppm_create_collection(DTYPE(part_mapping)_,DTYPE(part_mapping)_,generate="abstract")

#if __MYTYPE == __MAPPINGTYPE
      TYPE,ABSTRACT,EXTENDS(ppm_t_mapping_) :: ppm_t_mesh_mapping_
          !------------------------------------------------------------------
          !  Mesh ghosts mappings
          !------------------------------------------------------------------
          INTEGER,               DIMENSION(:),   POINTER :: ghostsize         => NULL()
          !!! size of the ghost layer, in number of mesh points
          INTEGER,               DIMENSION(:),   POINTER :: ghost_fromsub     => NULL()
          !!! list of source subs of ghost mesh blocks (globel sub number).
          !!! These are the owner subs of the actual real mesh points
          !!! 1st index: meshblock ID
          INTEGER,               DIMENSION(:),   POINTER :: ghost_tosub       => NULL()
          !!! list of target subs of ghost mesh blocks (globel sub number).
          !!! These are the subs a block will serve as a ghost on.
          !!! 1st index: meshblock ID
          !Yaser
          INTEGER,               DIMENSION(:,:), POINTER :: ghost_patchid     => NULL()
          !!! list of patches of ghost mesh blocks (globel sub number).
          !!! 1st index: patch ID
          INTEGER,               DIMENSION(:,:), POINTER :: ghost_blkstart    => NULL()
          !!! start (lower-left corner) of ghost mesh block in GLOBAL
          !!! mesh coordinates. First index: x,y[,z], 2nd: meshblock ID
          INTEGER,               DIMENSION(:,:), POINTER :: ghost_blksize     => NULL()
          !!! size (in grid points) of ghost blocks. 1st index: x,y[,z], 2nd:
          !!! meshblock ID
          INTEGER,               DIMENSION(:),   POINTER :: ghost_blk         => NULL()
          !!! mesh ghost block list. 1st index: target processor
          INTEGER                                        :: ghost_nsend
          !!! number of mesh blocks to be sent as ghosts
          INTEGER                                        :: ghost_nrecv
          !!! number of mesh blocks to be recvd as ghosts
          INTEGER,               DIMENSION(:),   POINTER :: ghost_recvtosub   => NULL()
          !!! list of target subs for ghost mesh blocks to be received,
          !!! i.e. being ghost on the local processor (globel sub number).
          !!! These are the subs where the blocks will serve as ghosts
          !!! 1st index: meshblock ID
          INTEGER,               DIMENSION(:,:), POINTER :: ghost_recvpatchid => NULL()
          !!! list of patches (global indices) for ghost mesh blocks to be received,
          !!! i.e. being ghost on the local processor (globel sub number).
          !!! 1st index: patch ID
          INTEGER,               DIMENSION(:,:), POINTER :: ghost_recvblkstart=> NULL()
          !!! start (lower-left corner) of received ghost mesh block in
          !!! GLOBAL  mesh coordinates. 1st index: x,y[,z], 2nd: meshblock ID
          INTEGER,               DIMENSION(:,:), POINTER :: ghost_recvblksize => NULL()
          !!! size (in grid points) of recvd ghost blocks.
          !!! 1st index: x,y[,z], 2nd: meshblock ID
          INTEGER,               DIMENSION(:),   POINTER :: ghost_recvblk     => NULL()
          !!! mesh ghost block receive list. 1st index: target processor

          !----------------------------------------------------------------------
          !  Mesh mapping, send and receive lists
          !----------------------------------------------------------------------
          INTEGER, DIMENSION(:  ), POINTER :: ppm_mesh_isendfromsub  => NULL()
          !!! list of source subs to send from local processor (local sub number
          !!! on source processor)
          INTEGER, DIMENSION(:,:), POINTER :: ppm_mesh_isendblkstart => NULL()
          ! start (lower-left corner) of mesh block to be sent in GLOBAL
          ! mesh coordinates. First index: x,y[,z], 2nd: isendlist
          INTEGER, DIMENSION(:,:), POINTER :: ppm_mesh_isendpatchid  => NULL()
          !!! list of source patch ids to send from local processor
          INTEGER, DIMENSION(:,:), POINTER :: ppm_mesh_isendblksize  => NULL()
          ! size (in grid points) of blocks to be sent

          INTEGER, DIMENSION(:  ), POINTER :: ppm_mesh_irecvtosub    => NULL()
          ! list of destination subs to recv to on local processors (local sub
          ! number on destination processor)
          INTEGER, DIMENSION(:,:), POINTER :: ppm_mesh_irecvblkstart => NULL()
          ! start (lower-left corner) of mesh block to be recvd in GLOBAL
          ! mesh coordinates. First index: x,y[,z], 2nd: isendlist
          INTEGER, DIMENSION(:,:), POINTER :: ppm_mesh_irecvpatchid  => NULL()
          !!! list of source patch ids to receive from local processor
          INTEGER, DIMENSION(:,:), POINTER :: ppm_mesh_irecvblksize  => NULL()
          ! size (in grid points) of blocks to be recvd

      CONTAINS
           PROCEDURE(map_mesh_create_),  DEFERRED :: create
           PROCEDURE(map_mesh_destroy_), DEFERRED :: destroy
!              PROCEDURE :: mesh_send
!             GENERIC :: mesh_pop => &
!             &          DTYPE(ppm_map_mesh_pop_1d),  &
!             &          DTYPE(ppm_map_mesh_pop_1dc), &
!             &          ppm_map_mesh_pop_1di,        &
!             &          ppm_map_mesh_pop_1dl,        &
!             &          DTYPE(ppm_map_mesh_pop_2d),  &
!             &          DTYPE(ppm_map_mesh_pop_2dc), &
!             &          ppm_map_mesh_pop_2di,        &
!             &          ppm_map_mesh_pop_2dl
!             GENERIC :: mesh_push =>                 &
!             &          DTYPE(ppm_map_mesh_push_1d), &
!             &          DTYPE(ppm_map_mesh_push_1dc),&
!             &          ppm_map_mesh_push_1di,       &
!             &          ppm_map_mesh_push_1dl,       &
!             &          DTYPE(ppm_map_mesh_push_2d), &
!             &          DTYPE(ppm_map_mesh_push_2dc),&
!             &          ppm_map_mesh_push_2di,       &
!             &          ppm_map_mesh_push_2dl
      END TYPE ppm_t_mesh_mapping_
minclude ppm_create_collection(mesh_mapping_,mesh_mapping_,generate="abstract")
#endif
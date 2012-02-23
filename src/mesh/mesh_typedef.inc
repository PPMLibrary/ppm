         TYPE DTYPE(ppm_t_mesh_maplist)
         !!! TODO: check what this is used for (imported from Petros code
             INTEGER, POINTER  :: target_topoid => NULL()
             !!! target topology ID
             INTEGER, POINTER  :: target_meshid => NULL()
             !!! target mesh ID
             INTEGER, POINTER  :: nsendlist => NULL()
             !!! send rank list size
             INTEGER, POINTER  :: nrecvlist => NULL()
             !!! recv rank size
             INTEGER, DIMENSION(:), POINTER  :: isendlist => NULL()
             !!! send rank lists
             INTEGER, DIMENSION(:), POINTER  :: irecvlist => NULL()
             !!! recv rank lists
             INTEGER, DIMENSION(:), POINTER  :: isendfromsub => NULL()
             !!! source send sublist
             INTEGER, DIMENSION(:,:), POINTER  :: isendblkstart => NULL()
             !!! send block start list
             INTEGER, DIMENSION(:,:), POINTER  :: isendblksize  => NULL()
             !!! send block size list
             INTEGER, DIMENSION(:), POINTER  :: irecvtosub   => NULL()
             !!! recv send sublist
             INTEGER, DIMENSION(:,:), POINTER  :: irecvblkstart => NULL()
             !!! recv block start list
             INTEGER, DIMENSION(:,:), POINTER  :: irecvblksize => NULL()
             !!! recv block size list
             INTEGER, DIMENSION(:), POINTER  :: psendbuffer => NULL()
             !!! send buffer pointer
             INTEGER, DIMENSION(:), POINTER  :: precvbuffer => NULL()
             !!! recv buffer pointer
         END TYPE DTYPE(ppm_t_mesh_maplist)

         TYPE DTYPE(ppm_t_equi_mesh)
         !!! Type for equispaced cartesian meshes on subs
             INTEGER                           :: ID
             !!! ID of the mesh in the belonging topology
             !!! It is the same as its index in the ppm_t_topo%mesh array
             INTEGER, DIMENSION(:,:), POINTER  :: nnodes => NULL()
             !!! The number of mesh *nodes* (not cells) in each direction in
             !!! each sub

             INTEGER, DIMENSION(:,:), POINTER  :: istart => NULL()
             !!! Starting indices of the mesh of this sub in the global mesh

             INTEGER, DIMENSION(:  ), POINTER  :: Nm    => NULL()
             !!! global number of mesh points in computational domain

             !------------------------------------------------------------------
             !  Mesh ghosts mappings
             !------------------------------------------------------------------
             LOGICAL                          :: ghost_initialized = .FALSE.
             !!! is .TRUE. if the ghost mappings have been initialized
             !!! else, .FALSE.
             INTEGER, DIMENSION(:),   POINTER :: ghost_fromsub => NULL()
             !!! list of source subs of ghost mesh blocks (globel sub number).
             !!! These are the owner subs of the actual real mesh points
             !!! 1st index: meshblock ID
             INTEGER, DIMENSION(:),   POINTER :: ghost_tosub   => NULL()
             !!! list of target subs of ghost mesh blocks (globel sub number).
             !!! These are the subs a block will serve as a ghost on.
             !!! 1st index: meshblock ID
             INTEGER, DIMENSION(:,:), POINTER :: ghost_blkstart => NULL()
             !!! start (lower-left corner) of ghost mesh block in GLOBAL
             !!! mesh coordinates. First index: x,y[,z], 2nd: meshblock ID
             INTEGER, DIMENSION(:,:), POINTER :: ghost_blksize  => NULL()
             !!! size (in grid points) of ghost blocks. 1st index: x,y[,z], 2nd:
             !!! meshblock ID
             INTEGER, DIMENSION(:)  , POINTER :: ghost_blk      => NULL()
             !!! mesh ghost block list. 1st index: target processor
             INTEGER                          :: ghost_nsend
             !!! number of mesh blocks to be sent as ghosts
             INTEGER                          :: ghost_nrecv
             !!! number of mesh blocks to be recvd as ghosts
             INTEGER, DIMENSION(:), POINTER   :: ghost_recvtosub => NULL()
             !!! list of target subs for ghost mesh blocks to be received,
             !!! i.e. being ghost on the local processor (globel sub number).
             !!! These are the subs where the blocks will serve as ghosts
             !!! 1st index: meshblock ID
             INTEGER, DIMENSION(:,:), POINTER  :: ghost_recvblkstart => NULL()
             !!! start (lower-left corner) of received ghost mesh block in
             !!! GLOBAL  mesh coordinates. 1st index: x,y[,z], 2nd: meshblock ID
             INTEGER, DIMENSION(:,:), POINTER  :: ghost_recvblksize => NULL()
             !!! size (in grid points) of recvd ghost blocks.
             !!! 1st index: x,y[,z], 2nd: meshblock ID
             INTEGER, DIMENSION(:)  , POINTER  :: ghost_recvblk => NULL()
             !!! mesh ghost block receive list. 1st index: target processor

             TYPE(DTYPE(ppm_t_mesh_maplist)), POINTER :: mapping => NULL()


         END TYPE DTYPE(ppm_t_equi_mesh)

#undef DTYPE
#undef MK

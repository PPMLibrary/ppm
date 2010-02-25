      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_data_mesh
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all data structures and
      !                 definitions that are PRIVATE to the mesh routines.
      !                 It also included those routines and provides
      !                 INTERFACEs.
      !                
      !  Remarks      : The terminology distinguishes between meshes and
      !                 fields (the data living on the meshes). Several
      !                 fields can use the same mesh. Meshes are defined as
      !                 ppm-internal TYPES, whereas fields are
      !                 user-provided arrays.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_mesh.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:28:18  ivos
      !  Initial implementation. Originated from splitting the old ppm
      !  modules.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_data_mesh

         !----------------------------------------------------------------------
         !  Mesh TYPEs
         !----------------------------------------------------------------------
         ! ppm-internal structure for equispaced cartesian meshes on subs
         TYPE ppm_type_equi_mesh
             ! topology ID (internal numbering) this mesh lives on
             INTEGER                           :: topoid
             ! The number of mesh NODES (not cells) in each direction in
             ! each sub
             INTEGER, DIMENSION(:,:), POINTER  :: nnodes
             ! Starting indices of the mesh of this sub in the global mesh
             INTEGER, DIMENSION(:,:), POINTER  :: istart
             ! global number of mesh points in computational domain
             INTEGER, DIMENSION(:  ), POINTER  :: Nm
         END TYPE

         ! mesh lists with translation between user numbering and internal
         ! numbering of mesh IDs
         TYPE ppm_type_mesh_list
             ! translation from internal to user numbers
             INTEGER, DIMENSION(:  ), POINTER  :: user
             ! translation from user numbers to internal ones
             INTEGER, DIMENSION(:  ), POINTER  :: internal
         END TYPE

         !----------------------------------------------------------------------
         !  Internal meshes
         !----------------------------------------------------------------------
         ! first index: meshid (idea: several meshes can be defined on the same
         ! topology at once, e.g. for multigrid solvers.), second: topoid
         ! (internal numbering)
         TYPE(ppm_type_equi_mesh), DIMENSION(:,:), POINTER :: ppm_cart_mesh
         ! highest internal mesh ID available on each topology (internal
         ! numbering)
         INTEGER, DIMENSION(:  ), POINTER                  :: ppm_max_meshid
         ! number translation lists (user<-->internal) for each topology
         TYPE(ppm_type_mesh_list), DIMENSION(:  ), POINTER :: ppm_meshid

         !----------------------------------------------------------------------
         !  Mesh mapping, send and receive lists
         !----------------------------------------------------------------------
         ! list of source subs to send from local processor (local sub number
         ! on source processor)
         INTEGER, DIMENSION(:  ), POINTER          :: ppm_mesh_isendfromsub
         ! start (lower-left corner) of mesh block to be sent in GLOBAL
         ! mesh coordinates. First index: x,y[,z], 2nd: isendlist
         INTEGER, DIMENSION(:,:), POINTER          :: ppm_mesh_isendblkstart
         ! size (in grid points) of blocks to be sent
         INTEGER, DIMENSION(:,:), POINTER          :: ppm_mesh_isendblksize
         ! list of destination subs to recv to on local processors (local sub
         ! number on destination processor)
         INTEGER, DIMENSION(:  ), POINTER          :: ppm_mesh_irecvtosub
         ! start (lower-left corner) of mesh block to be recvd in GLOBAL
         ! mesh coordinates. First index: x,y[,z], 2nd: isendlist
         INTEGER, DIMENSION(:,:), POINTER          :: ppm_mesh_irecvblkstart
         ! size (in grid points) of blocks to be recvd
         INTEGER, DIMENSION(:,:), POINTER          :: ppm_mesh_irecvblksize

         ! define the target mesh id for mapping
         INTEGER          :: ppm_target_meshid = -1
         ! define the source mesh id for mapping
         INTEGER          :: ppm_source_meshid = -1
         ! define the target mesh id for the mapping 
         ! NOTICE: NOT the same variable as the ppm_target_topoid in the
         ! module ppm_module_map
         INTEGER          :: ppm_target_topoid = -1

         !----------------------------------------------------------------------
         !  Mesh ghosts mappings
         !----------------------------------------------------------------------
         ! list of source subs of ghost mesh blocks (globel sub number).
         ! These are the owner subs of the actual real mesh points
         ! 1st index: meshblock ID, 2nd: meshid (internal numbering), 3rd:
         ! topoid (internal numbering)
         INTEGER, DIMENSION(:,:,:), POINTER            :: ppm_mesh_ghost_fromsub
         ! list of target subs of ghost mesh blocks (globel sub number).
         ! These are the subs a block will serve as a ghost on.
         ! 1st index: meshblock ID, 2nd: meshid (internal numbering), 3rd:
         ! topoid (internal numbering)
         INTEGER, DIMENSION(:,:,:), POINTER            :: ppm_mesh_ghost_tosub
         ! start (lower-left corner) of ghost mesh block in GLOBAL
         ! mesh coordinates. First index: x,y[,z], 2nd: meshblock ID, 3rd:
         ! meshid (internal numbering), 4th: topoid (internal numbering)
         INTEGER, DIMENSION(:,:,:,:), POINTER          ::    &
     &       ppm_mesh_ghost_blkstart
         ! size (in grid points) of ghost blocks. 1st index: x,y[,z], 2nd:
         ! meshblock ID, 3rd: meshid (internal numbering), 4th: topoid
         ! (internal numbering)
         INTEGER, DIMENSION(:,:,:,:), POINTER          :: ppm_mesh_ghost_blksize
         ! mesh ghost block list. 1st index: target processor, 2nd: meshid
         ! (internal numbering), 3rd: topoid (internal numbering)
         INTEGER, DIMENSION(:,:,:)  , POINTER          :: ppm_mesh_ghost_blk
         ! number of mesh blocks to be sent as ghosts. 1st: meshid, 2nd:
         ! topoid (both internal numbering)
         INTEGER, DIMENSION(:,:)    , POINTER          :: ppm_mesh_ghost_nsend
         ! number of mesh blocks to be recvd as ghosts. 1st: meshid, 2nd:
         ! topoid (both internal numbering)
         INTEGER, DIMENSION(:,:)    , POINTER          :: ppm_mesh_ghost_nrecv
         ! list of target subs for ghost mesh blocks to be received,
         ! i.e. being ghost on the local processor (globel sub number).
         ! These are the subs where the blocks will serve as ghosts
         ! 1st index: meshblock ID, 2nd: meshid (internal numbering), 3rd:
         ! topoid (internal numbering)
         INTEGER, DIMENSION(:,:,:), POINTER            ::    &
     &       ppm_mesh_ghost_recvtosub
         ! start (lower-left corner) of received ghost mesh block in GLOBAL
         ! mesh coordinates. First index: x,y[,z], 2nd: meshblock ID, 3rd:
         ! meshid (internal numbering), 4th: topoid (internal numbering)
         INTEGER, DIMENSION(:,:,:,:), POINTER          ::    &
     &       ppm_mesh_ghost_recvblkstart
         ! size (in grid points) of recvd ghost blocks. 1st index: x,y[,z], 2nd:
         ! meshblock ID, 3rd: meshid (internal numbering), 4th: topoid
         ! (internal numbering)
         INTEGER, DIMENSION(:,:,:,:), POINTER          ::    &
     &       ppm_mesh_ghost_recvblksize
         ! mesh ghost block receive list. 1st index: target processor, 2nd: 
         ! meshid (internal numbering), 3rd: topoid (internal numbering)
         INTEGER, DIMENSION(:,:,:)  , POINTER          :: ppm_mesh_ghost_recvblk

      END MODULE ppm_module_data_mesh

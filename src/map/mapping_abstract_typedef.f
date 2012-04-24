!----------------------------------------------------------------------
!  Generic mapping data type
!----------------------------------------------------------------------
TYPE,ABSTRACT :: DTYPE(ppm_t_mapping)_
    INTEGER :: source_topoid
    !!! topology ID on which the data was stored before computing the
    !!! mapping 
    INTEGER :: target_topoid
    !!! topology ID on which the data should be stored after the
    !!! mapping

    !------------------------------------------------------------------
    !  buffers for communication
    !------------------------------------------------------------------
    REAL(MK), DIMENSION(:), POINTER  :: send => NULL()
    REAL(MK), DIMENSION(:), POINTER  :: recv => NULL()
    INTEGER, DIMENSION(:),  POINTER  :: nsend => NULL()
    INTEGER, DIMENSION(:),  POINTER  :: nrecv => NULL()
    INTEGER, DIMENSION(:),  POINTER  :: psend => NULL()
    INTEGER, DIMENSION(:),  POINTER  :: precv => NULL()
    INTEGER, DIMENSION(:,:),POINTER  :: pp    => NULL()
    INTEGER, DIMENSION(:,:),POINTER  :: qq    => NULL()
    INTEGER                          :: old_nsendlist = 0
    INTEGER                          :: old_buffer_set = 0

    REAL(MK),DIMENSION(:),POINTER :: ppm_recvbuffer => NULL()
    !!! receive buffer for particles 

    REAL(MK),DIMENSION(:),POINTER :: ppm_sendbuffer => NULL()
    !!! send buffer for particles 

    INTEGER            ,DIMENSION(:,:),POINTER :: ppm_ghosthack => NULL()
    !!! invert map of ghost for symmetry

    INTEGER             ,DIMENSION(:),POINTER  :: ppm_psendbuffer => NULL()
    !!! pointer to particles within the send buffer
    !!!
    !!! In terms of particle *not* the actual position in the buffer
    INTEGER             ,DIMENSION(:),POINTER  :: ppm_precvbuffer => NULL()
    !!! pointer to particles within the recv buffer
    !!!
    !!! In terms of particle *not* the actual position in the buffer

    INTEGER                                        ::  ppm_nsendbuffer
    !!! the size of the send buffer (not the array but the number of
    !!! used elements)
    !!!
    !!! In terms of entries in the buffer *not* the number of particles 
    INTEGER                                        ::  ppm_nrecvbuffer
    !!! the size of the recv buffer (not the array but the number of
    !!! used elements)
    !!!
    !!! In terms of entries in the buffer *not* the number of particles 
    INTEGER                                        ::  ppm_sendbufsize
    !!! the actual size of the send buffer array
    INTEGER                                        ::  ppm_recvbufsize
    !!! the actual size of the receive buffer array
    INTEGER                                        ::  ppm_buffer_set
    !!! the total number of particle fields packed in
    !!! the send buffer, ie. xp, vp is two sets

    INTEGER            , DIMENSION(:),POINTER :: ppm_buffer2part => NULL()
    !!! Used for the original on-processor particle IDs in the order in which
    !!! they are in the sendbuffer. Used to push additional particle
    !!! data on the buffer in the correct order.
    INTEGER            , DIMENSION(:),POINTER :: ppm_buffer_type => NULL()
    !!! types of the data on the sendbuffer
    INTEGER            , DIMENSION(:),POINTER :: ppm_buffer_dim => NULL()
    !!! dimensions of the original on-processor particle arrays. 

    REAL(MK),DIMENSION(:),POINTER             :: ppm_ghost_offset => NULL()
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

    !----------------------------------------------------------------------
    !  mapping
    !----------------------------------------------------------------------
    INTEGER                        :: ppm_map_type 
    INTEGER                        :: ppm_nsendlist
    INTEGER                        :: ppm_nrecvlist
    INTEGER, DIMENSION(:), POINTER :: ppm_isendlist => NULL()
    INTEGER, DIMENSION(:), POINTER :: ppm_irecvlist => NULL()
END TYPE DTYPE(ppm_t_mapping)_

TYPE,ABSTRACT, EXTENDS(DTYPE(ppm_t_mapping)_) :: DTYPE(ppm_t_part_mapping)_
    INTEGER :: oldNpart
    INTEGER :: newNpart

    CONTAINS

    PROCEDURE(DTYPE(map_create)_),DEFERRED     :: create
    PROCEDURE(DTYPE(map_destroy)_),DEFERRED     :: destroy

END TYPE DTYPE(ppm_t_part_mapping)_
minclude define_abstract_collection_type(DTYPE(ppm_t_part_mapping)_)


TYPE,ABSTRACT,EXTENDS(DTYPE(ppm_t_mapping)_) :: DTYPE(ppm_t_mesh_mapping)_

    !             
    !             PROCEDURE       :: mesh_send
    !             GENERIC, PUBLIC ::  mesh_pop => &
    !                    DTYPE(ppm_map_mesh_pop_1d),&
    !                    DTYPE(ppm_map_mesh_pop_1dc),&
    !                    ppm_map_mesh_pop_1di,&
    !                    ppm_map_mesh_pop_1dl,&
    !                    DTYPE(ppm_map_mesh_pop_2d),&
    !                    DTYPE(ppm_map_mesh_pop_2dc),&
    !                    ppm_map_mesh_pop_2di,&
    !                    ppm_map_mesh_pop_2dl
    !             GENERIC, PUBLIC ::  mesh_push => &
    !                    DTYPE(ppm_map_mesh_push_1d),&
    !                    DTYPE(ppm_map_mesh_push_1dc),&
    !                    ppm_map_mesh_push_1di,&
    !                    ppm_map_mesh_push_1dl,&
    !                    DTYPE(ppm_map_mesh_push_2d),&
    !                    DTYPE(ppm_map_mesh_push_2dc),&
    !                    ppm_map_mesh_push_2di,&
    !                    ppm_map_mesh_push_2dl
END TYPE DTYPE(ppm_t_mesh_mapping)_

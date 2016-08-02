#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2
#define __INTEGER                  3
#define __LOGICAL                  4
#define __SINGLE_PRECISION_COMPLEX 5
#define __DOUBLE_PRECISION_COMPLEX 6
#define __MAPPINGTYPE              9

      MODULE ppm_module_mapping_typedef
      !!! This module defines the mapping data type and provides
      !!!  the basic mapping routines for particles and meshes; namely
      !!!  the push, send and pop routine.
      !----------------------------------------------------------------------
      !  Modules
      !----------------------------------------------------------------------
      !   USE ppm_module_map_part_util
      !   USE ppm_module_map_part_ghost
      !   USE ppm_module_map_part_global
      !   USE ppm_module_map_part_partial
      USE ppm_module_data, ONLY : ppm_kind_double,ppm_kind_single,ppm_char, &
      &   ppm_error_error,ppm_kind_int64
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_interfaces, ONLY : ppm_t_mapping_,ppm_t_part_mapping_s_, &
      &   ppm_t_part_mapping_d_,ppm_c_part_mapping_s_,ppm_c_part_mapping_d_,  &
      &   ppm_t_mesh_mapping_,ppm_c_mesh_mapping_,ppm_t_ptr_part_mapping_s,   &
      &   ppm_t_ptr_part_mapping_d,ppm_t_ptr_mesh_mapping
      IMPLICIT NONE

      PRIVATE

      !----------------------------------------------------------------------
      !  buffers for communication
      !----------------------------------------------------------------------
      REAL(ppm_kind_single), DIMENSION(:),   POINTER :: ppm_sendbuffers => NULL()
      REAL(ppm_kind_double), DIMENSION(:),   POINTER :: ppm_sendbufferd => NULL()
      !!! send buffer
      REAL(ppm_kind_single), DIMENSION(:),   POINTER :: ppm_recvbuffers => NULL()
      REAL(ppm_kind_double), DIMENSION(:),   POINTER :: ppm_recvbufferd => NULL()
      !!! recv buffer


      INTEGER,               DIMENSION(:,:), POINTER :: ppm_ghosthack   => NULL()
      !!! invert map of ghost for symmetry

      INTEGER,               DIMENSION(:),   POINTER :: ppm_psendbuffer => NULL()
      !!! pointer to particles within the send buffer
      !!!
      !!! In terms of particle *not* the actual position in the buffer
      INTEGER,               DIMENSION(:),   POINTER :: ppm_precvbuffer => NULL()
      !!! pointer to particles within the recv buffer
      !!!
      !!! In terms of particle *not* the actual position in the buffer
      INTEGER(ppm_kind_int64)                        :: ppm_nsendbuffer
      !!! the size of the send buffer (not the array but the number of
      !!! used elements)
      !!!
      !!! In terms of entries in the buffer *not* the number of particles
      INTEGER(ppm_kind_int64)                        :: ppm_nrecvbuffer
      !!! the size of the recv buffer (not the array but the number of
      !!! used elements)
      !!!
      !!! In terms of entries in the buffer *not* the number of particles
      INTEGER,               DIMENSION(:),   POINTER :: ppm_buffer2part => NULL()
      !!! Used for the original on-processor particle IDs in the order in which
      !!! they are in the sendbuffer. Used to push additional particle
      !!! data on the buffer in the correct order.
      INTEGER                                        :: ppm_buffer_set
      !!! the total number of particle fields packed in
      !!! the send buffer, ie. xp, vp is two sets
      INTEGER,               DIMENSION(:),   POINTER :: ppm_buffer_type => NULL()
      !!! types of the data on the sendbuffer
      INTEGER,               DIMENSION(:),   POINTER :: ppm_buffer_dim => NULL()
      !!! dimensions of the original on-processor particle arrays.
      REAL(ppm_kind_single), DIMENSION(:),   POINTER :: ppm_ghost_offsets => NULL()
      REAL(ppm_kind_double), DIMENSION(:),   POINTER :: ppm_ghost_offsetd => NULL()
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
      !!! `ppm_ghost_offset(ibuffer+0) = xp_offset(1,)`               +
      !!! `ppm_ghost_offset(ibuffer+1) = xp_offset(2,)`               +
      !!! `ppm_ghost_offset(ibuffer+2) = xp_offset(3,)`
      REAL(ppm_kind_single), DIMENSION(:),   POINTER :: ppm_ghost_offset_facs => NULL()
      REAL(ppm_kind_double), DIMENSION(:),   POINTER :: ppm_ghost_offset_facd => NULL()
      !!! ghost offset factor
      !!!
      !!! This buffer is needed for the special case of symmetric and
      !!! antisymmetric boundary conditions, because the offset is not a
      !!! simple constant to be added to the particle coordinate

      !----------------------------------------------------------------------
      !  mapping
      !----------------------------------------------------------------------
      INTEGER                                        :: ppm_map_type
      INTEGER                                        :: ppm_nsendlist
      INTEGER                                        :: ppm_nrecvlist
      INTEGER,               DIMENSION(:),   POINTER :: ppm_isendlist => NULL()
      INTEGER,               DIMENSION(:),   POINTER :: ppm_irecvlist => NULL()

      LOGICAL,               DIMENSION(:),   POINTER :: ppm_lsendlist => NULL()
      LOGICAL,               DIMENSION(:),   POINTER :: ppm_lrecvlist => NULL()
      !----------------------------------------------------------------------
      !!! Holds mesh mapping, receive buffers and
      !!! mesh ghost layer mapping buffers.
      !!!
      !!! [NOTE]
      !!! These variables declared in should not be accessed by the
      !!! PPM client developer. They are managed interally by the library.
      !----------------------------------------------------------------------
      !  Mesh mapping, send and receive lists
      !----------------------------------------------------------------------
      INTEGER,               DIMENSION(:  ), POINTER :: ppm_mesh_isendfromsub  => NULL()
      !!! list of source subs to send from local processor (local sub number
      !!! on source processor)
      INTEGER,               DIMENSION(:,:), POINTER :: ppm_mesh_isendblkstart => NULL()
      ! start (lower-left corner) of mesh block to be sent in GLOBAL
      ! mesh coordinates. First index: x,y[,z], 2nd: isendlist
      INTEGER,               DIMENSION(:,:), POINTER :: ppm_mesh_isendpatchid  => NULL()
      !!! list of source patch ids to send from local processor
      INTEGER,               DIMENSION(:,:), POINTER :: ppm_mesh_isendblksize  => NULL()
      ! size (in grid points) of blocks to be sent

      INTEGER,               DIMENSION(:  ), POINTER :: ppm_mesh_irecvtosub    => NULL()
      ! list of destination subs to recv to on local processors (local sub
      ! number on destination processor)
      INTEGER,               DIMENSION(:,:), POINTER :: ppm_mesh_irecvblkstart => NULL()
      ! start (lower-left corner) of mesh block to be recvd in GLOBAL
      ! mesh coordinates. First index: x,y[,z], 2nd: isendlist
      INTEGER,               DIMENSION(:,:), POINTER :: ppm_mesh_irecvpatchid  => NULL()
      !!! list of source patch ids to receive from local processor
      INTEGER,               DIMENSION(:,:), POINTER :: ppm_mesh_irecvblksize  => NULL()
      ! size (in grid points) of blocks to be recvd

#ifdef __MPI
      INTEGER,               DIMENSION(:), ALLOCATABLE :: ppm_request
      INTEGER                                          :: ppm_nrequest
#endif

      PUBLIC :: ppm_sendbuffers
      PUBLIC :: ppm_sendbufferd
      PUBLIC :: ppm_recvbuffers
      PUBLIC :: ppm_recvbufferd

      PUBLIC :: ppm_ghosthack
      PUBLIC :: ppm_psendbuffer
      PUBLIC :: ppm_precvbuffer
      PUBLIC :: ppm_nsendbuffer
      PUBLIC :: ppm_nrecvbuffer
      PUBLIC :: ppm_buffer_set
      PUBLIC :: ppm_buffer2part
      PUBLIC :: ppm_buffer_type
      PUBLIC :: ppm_buffer_dim
      PUBLIC :: ppm_ghost_offsets
      PUBLIC :: ppm_ghost_offsetd
      PUBLIC :: ppm_ghost_offset_facs
      PUBLIC :: ppm_ghost_offset_facd
      PUBLIC :: ppm_map_type
      PUBLIC :: ppm_nsendlist
      PUBLIC :: ppm_nrecvlist
      PUBLIC :: ppm_isendlist
      PUBLIC :: ppm_irecvlist
      PUBLIC :: ppm_lsendlist
      PUBLIC :: ppm_lrecvlist

      PUBLIC :: ppm_mesh_isendfromsub
      PUBLIC :: ppm_mesh_isendblkstart
      PUBLIC :: ppm_mesh_isendpatchid
      PUBLIC :: ppm_mesh_isendblksize
      PUBLIC :: ppm_mesh_irecvtosub
      PUBLIC :: ppm_mesh_irecvblkstart
      PUBLIC :: ppm_mesh_irecvpatchid
      PUBLIC :: ppm_mesh_irecvblksize

#ifdef __MPI
      PUBLIC :: ppm_request
      PUBLIC :: ppm_nrequest
#endif
      !----------------------------------------------------------------------
      !  Type declaration
      !----------------------------------------------------------------------
#define  DTYPE(a) a/**/_s
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#define __MYTYPE __MAPPINGTYPE
#include "map/mapping_typedef.f"

#define  DTYPE(a) a/**/_d
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "map/mapping_typedef.f"

      PUBLIC :: ppm_t_mapping

      PUBLIC :: ppm_t_part_mapping_s
      PUBLIC :: ppm_t_part_mapping_d
      PUBLIC :: ppm_c_part_mapping_s
      PUBLIC :: ppm_c_part_mapping_d

      PUBLIC :: ppm_t_mesh_mapping
      PUBLIC :: ppm_c_mesh_mapping

      !----------------------------------------------------------------------
      !  Type-bound procedures
      !----------------------------------------------------------------------
      CONTAINS

#define  DTYPE(a) a/**/_s
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#define __MYTYPE __MAPPINGTYPE
#include "map/mapping_typeproc.f"

#define  DTYPE(a) a/**/_d
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "map/mapping_typeproc.f"

      END MODULE ppm_module_mapping_typedef

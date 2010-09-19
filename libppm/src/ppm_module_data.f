      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_data
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Lab (ETH Zurich), 
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

      MODULE ppm_module_data
      !!! Declares global data types and variables.
      !!!
      !!! [NOTE]
      !!! Most of the declared variables in this module should not be accessed
      !!! directly by the PPM client developer, they are used internally in the
      !!! library.
         
         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_typedef

         !----------------------------------------------------------------------
         !  Global TYPEs
         !----------------------------------------------------------------------

         !----------------------------------------------------------------------
         !  buffers for communication
         !----------------------------------------------------------------------
         REAL(ppm_kind_double), DIMENSION(  :), POINTER ::  ppm_sendbufferd 
         !!! send buffer for particles 
         REAL(ppm_kind_double), DIMENSION(  :), POINTER ::  ppm_recvbufferd
         !!! recv buffer for particles

         REAL(ppm_kind_single), DIMENSION(  :), POINTER ::  ppm_sendbuffers
         !!! send buffer for particles
         REAL(ppm_kind_single), DIMENSION(  :), POINTER ::  ppm_recvbuffers
         !!! recv buffer for particles

         INTEGER              , DIMENSION(:,:), POINTER ::  ppm_ghosthack
         !!! invert map of ghost for symmetry

         INTEGER              , DIMENSION(  :), POINTER ::  ppm_psendbuffer
         !!! pointer to particles within the send buffer
         !!!
         !!! In terms of particle *not* the actual position in the buffer
         INTEGER              , DIMENSION(  :), POINTER ::  ppm_precvbuffer
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

         INTEGER              , DIMENSION(  :), POINTER ::  ppm_buffer2part
         !!! Used for the original on-processor particle IDs in the order in which
         !!! they are in the sendbuffer. Used to push additional particle
         !!! data on the buffer in the correct order.
         INTEGER              , DIMENSION(  :), POINTER ::  ppm_buffer_type
         !!! types of the data on the sendbuffer
         INTEGER              , DIMENSION(  :), POINTER ::  ppm_buffer_dim
         !!! dimensions of the original on-processor particle arrays. 

         REAL(ppm_kind_single), DIMENSION(:), POINTER   ::  ppm_ghost_offsets
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
         REAL(ppm_kind_double), DIMENSION(:), POINTER  ::  ppm_ghost_offsetd
         !!! ghost offset (double precison)
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
         INTEGER, DIMENSION(:), POINTER :: ppm_isendlist
         INTEGER, DIMENSION(:), POINTER :: ppm_irecvlist

         !----------------------------------------------------------------------
         !  Precision
         !----------------------------------------------------------------------
         INTEGER :: ppm_kind
         !!! Precision
         INTEGER :: ppm_mpi_kind
         !!! MPI Precision

         !----------------------------------------------------------------------
         !  Dimensionality
         !----------------------------------------------------------------------
         INTEGER :: ppm_dim
         !!! Dimensionality

         !----------------------------------------------------------------------
         ! Topologies
         !----------------------------------------------------------------------
         TYPE(ppm_ptr_t_topo), DIMENSION(:), POINTER :: ppm_topo
         !!! the PPM topologies array

         INTEGER :: ppm_next_avail_topo
         !!! ID of the next available topology to be used by
         !!! ppm_topo_alloc.
         !!!
         !!! At initialization this is set to ppm_param_undefined              +
         !!! If it points within (1,SIZE(ppm_topo)) this slot is (re)used to
         !!! store the new topology, if it is > SIZE(ppm_topo) then the ppm_topo
         !!! array must be extended
         
         !----------------------------------------------------------------------
         !  Debugging
         !----------------------------------------------------------------------
         INTEGER :: ppm_debug = 0
         !!! are we in Debugging mode?

         !----------------------------------------------------------------------
         !  Has ppm_init been called?
         !----------------------------------------------------------------------
         LOGICAL :: ppm_initialized = .FALSE.
         !!! Has ppm_init been called?

         !----------------------------------------------------------------------
         !  parallel variables
         !----------------------------------------------------------------------
         INTEGER :: ppm_nproc
         INTEGER :: ppm_rank
         INTEGER :: ppm_comm
         ! relative speeds of the processors (for load balancing)
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: ppm_proc_speed

         !----------------------------------------------------------------------
         !  Numerical tolerance. Differences smaller than this are considered
         !  zero.
         !----------------------------------------------------------------------
         REAL(ppm_kind_double) :: ppm_myepsd
         REAL(ppm_kind_single) :: ppm_myepss

         !----------------------------------------------------------------------
         !  Constants (computed in ppm_init)
         !----------------------------------------------------------------------
         REAL(ppm_kind_double) :: ppm_pi_d
         REAL(ppm_kind_single) :: ppm_pi_s

         !----------------------------------------------------------------------
         !  I/O Units
         !----------------------------------------------------------------------
         INTEGER               :: ppm_stdout = 6
         INTEGER               :: ppm_stderr = 0
         INTEGER               :: ppm_logfile = -1
         
      END MODULE ppm_module_data

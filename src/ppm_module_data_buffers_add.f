      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_data_buffers_add
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

      MODULE ppm_module_data_buffers_add
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
         TYPE ppm_t_part_modify
             INTEGER                            :: Nrnew
             !number of new real particles added by the user
             INTEGER                            :: Ngnew
             !number of new ghost particles added by the user
             !INTEGER                            :: Ngsendnew
             !!number of new real particles that are ghosts for other procs
             !INTEGER                            :: Ngrecvnew
             !number of new ghost particles that are received from other procs
             INTEGER, DIMENSION(:), POINTER     :: idx_real_new => NULL()
             !indeces of new real particles 
             INTEGER, DIMENSION(:), POINTER     :: idx_ghost_new => NULL()
             !indeces of new ghost particles 
             !INTEGER, DIMENSION(:), POINTER     :: idx_ghost_send_new => NULL()
             !!indeces of new real particles that are ghosts for other procs
             INTEGER                             :: Npart
             !number of particles currently on this proc
             INTEGER                             :: Mpart
             !number of particles currently on this proc (incl. ghost parts)
             INTEGER                             :: Npart_add
             !number of particles being added on this proc
             INTEGER                             :: Npart_new
             !new total number of particles on this proc
             INTEGER                             :: Mpart_new
             !new total number of particles on this proc (incl. ghost parts)
         END TYPE

         !----------------------------------------------------------------------
         !  data
         !----------------------------------------------------------------------

         TYPE(ppm_t_part_modify), SAVE           :: modify


         !----------------------------------------------------------------------
         !  buffers for communication
         !----------------------------------------------------------------------
         REAL(ppm_kind_double),DIMENSION(:),POINTER :: ppm_sendbufferd => NULL()
         !!! send buffer for particles 
         REAL(ppm_kind_double),DIMENSION(:),POINTER :: ppm_recvbufferd => NULL()
         !!! recv buffer for particles

         REAL(ppm_kind_single),DIMENSION(:),POINTER :: ppm_sendbuffers => NULL()
         !!! send buffer for particles
         REAL(ppm_kind_single),DIMENSION(:),POINTER :: ppm_recvbuffers => NULL()
         !!! recv buffer for particles

         INTEGER            ,DIMENSION(:,:),POINTER :: ppm_ghosthack => NULL()
         !!! invert map of ghost for symmetry

         INTEGER             ,DIMENSION(:),POINTER :: ppm_psendbuffer => NULL()
         !!! pointer to particles within the send buffer
         !!!
         !!! In terms of particle *not* the actual position in the buffer
         INTEGER             ,DIMENSION(:),POINTER :: ppm_precvbuffer => NULL()
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

         REAL(ppm_kind_single),DIMENSION(:),POINTER::ppm_ghost_offsets => NULL()
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
         REAL(ppm_kind_double),DIMENSION(:),POINTER::ppm_ghost_offsetd => NULL()
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

         INTEGER                          :: old_nsendlist = 0
         INTEGER                          :: old_buffer_set = 0

         INTEGER                          :: add_mode = -1
         !!! true if it is REAL particles that are being popped out of 
         !!! the buffer. (the normal behaviour is false).

         INTEGER, PARAMETER               :: ppm_param_add_real_particles = 1
         INTEGER, PARAMETER               :: ppm_param_add_ghost_particles = 2


      END MODULE ppm_module_data_buffers_add

      !--*- f90 -*--------------------------------------------------------------
      !  Module   :                   ppm_module_typedef
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

      MODULE ppm_module_typedef
      !!! This module contains the definition of some ppm data types.
      !!! Namely the topology datatype.

         !----------------------------------------------------------------------
         !  Header file for global parameters
         !----------------------------------------------------------------------
         INCLUDE 'ppm_param.h'


         TYPE ppm_t_equi_mesh
         !!! Type for equispaced cartesian meshes on subs
             INTEGER                           :: ID
             !!! ID of the mesh in the belonging topology
             !!! It is the same as its index in the ppm_t_topo%mesh array
             INTEGER, DIMENSION(:,:), POINTER  :: nnodes
             !!! The number of mesh *nodes* (not cells) in each direction in
             !!! each sub

             INTEGER, DIMENSION(:,:), POINTER  :: istart
             !!! Starting indices of the mesh of this sub in the global mesh

             INTEGER, DIMENSION(:  ), POINTER  :: Nm
             !!! global number of mesh points in computational domain

             !------------------------------------------------------------------
             !  Mesh ghosts mappings
             !------------------------------------------------------------------
             LOGICAL                          :: ghost_initialized = .FALSE.
             !!! is .TRUE. if the ghost mappings have been initialized
             !!! else, .FALSE.
             INTEGER, DIMENSION(:),   POINTER :: ghost_fromsub
             !!! list of source subs of ghost mesh blocks (globel sub number).
             !!! These are the owner subs of the actual real mesh points
             !!! 1st index: meshblock ID
             INTEGER, DIMENSION(:),   POINTER :: ghost_tosub
             !!! list of target subs of ghost mesh blocks (globel sub number).
             !!! These are the subs a block will serve as a ghost on.
             !!! 1st index: meshblock ID
             INTEGER, DIMENSION(:,:), POINTER :: ghost_blkstart
             !!! start (lower-left corner) of ghost mesh block in GLOBAL
             !!! mesh coordinates. First index: x,y[,z], 2nd: meshblock ID
             INTEGER, DIMENSION(:,:), POINTER :: ghost_blksize
             !!! size (in grid points) of ghost blocks. 1st index: x,y[,z], 2nd:
             !!! meshblock ID
             INTEGER, DIMENSION(:)  , POINTER :: ghost_blk
             !!! mesh ghost block list. 1st index: target processor
             INTEGER                          :: ghost_nsend
             !!! number of mesh blocks to be sent as ghosts
             INTEGER                          :: ghost_nrecv
             !!! number of mesh blocks to be recvd as ghosts
             INTEGER, DIMENSION(:), POINTER   :: ghost_recvtosub
             !!! list of target subs for ghost mesh blocks to be received,
             !!! i.e. being ghost on the local processor (globel sub number).
             !!! These are the subs where the blocks will serve as ghosts
             !!! 1st index: meshblock ID
             INTEGER, DIMENSION(:,:), POINTER  :: ghost_recvblkstart
             !!! start (lower-left corner) of received ghost mesh block in
             !!! GLOBAL  mesh coordinates. 1st index: x,y[,z], 2nd: meshblock ID
             INTEGER, DIMENSION(:,:), POINTER  :: ghost_recvblksize
             !!! size (in grid points) of recvd ghost blocks.
             !!! 1st index: x,y[,z], 2nd: meshblock ID
             INTEGER, DIMENSION(:)  , POINTER  :: ghost_recvblk
             !!! mesh ghost block receive list. 1st index: target processor

             TYPE(ppm_t_mesh_maplist), POINTER :: mapping


         END TYPE

         TYPE ppm_t_mesh_maplist
         !!! TODO: check what this is used for (imported from Petros code
             INTEGER, POINTER  :: target_topoid
             !!! target topology ID
             INTEGER, POINTER  :: target_meshid
             !!! target mesh ID
             INTEGER, POINTER  :: nsendlist
             !!! send rank list size
             INTEGER, POINTER  :: nrecvlist
             !!! recv rank size
             INTEGER, DIMENSION(:), POINTER  :: isendlist
             !!! send rank lists
             INTEGER, DIMENSION(:), POINTER  :: irecvlist
             !!! recv rank lists
             INTEGER, DIMENSION(:), POINTER  :: isendfromsub
             !!! source send sublist
             INTEGER, DIMENSION(:,:), POINTER  :: isendblkstart
             !!! send block start list
             INTEGER, DIMENSION(:,:), POINTER  :: isendblksize
             !!! send block size list
             INTEGER, DIMENSION(:), POINTER  :: irecvtosub
             !!! recv send sublist
             INTEGER, DIMENSION(:,:), POINTER  :: irecvblkstart
             !!! recv block start list
             INTEGER, DIMENSION(:,:), POINTER  :: irecvblksize
             !!! recv block size list
             INTEGER, DIMENSION(:), POINTER  :: psendbuffer
             !!! send buffer pointer
             INTEGER, DIMENSION(:), POINTER  :: precvbuffer
             !!! recv buffer pointer
         END TYPE


         !----------------------------------------------------------------------
         !  Topology TYPE
         !----------------------------------------------------------------------
         TYPE ppm_t_topo
         !!! The topology type

            INTEGER                                      :: ID
            !!! ID of this topology
            !!!
            !!! It is the same as its index in the ppm_topo array
            LOGICAL                                      :: isdefined = .FALSE.
            !!! flag to tell if this topology is defined/in use
            INTEGER                                      :: prec
            !!! numerical precision (ppm_kind) for this topology

            REAL(ppm_kind_single), DIMENSION(:), POINTER :: min_physs
            !!! minimum of physical extend of the computational domain (single)
            !!! Note: first index is ppm_dim
            REAL(ppm_kind_single), DIMENSION(:), POINTER :: max_physs
            !!! maximum of physical extend of the computational domain (single)
            !!! Note: first index is ppm_dim
            REAL(ppm_kind_double), DIMENSION(:), POINTER :: min_physd
            !!! minimum of physical extend of the computational domain (double)
            !!! Note: first index is ppm_dim
            REAL(ppm_kind_double), DIMENSION(:), POINTER :: max_physd
            !!! maximum of physical extend of the computational domain (double)
            !!! Note: first index is ppm_dim

            INTEGER              , DIMENSION(:  ), POINTER :: bcdef
            !!! boundary conditions for the topology
            !!! Note: first index is 1-6 (each of the faces)

            INTEGER                                        :: nsubs
            !!! total number of subs on all processors.

            REAL(ppm_kind_single), DIMENSION(:,:), POINTER :: min_subs
            !!! mimimum of extension of subs (single)
            !!! Note: 1st index: x,y,(z), 2nd: subID
            REAL(ppm_kind_single), DIMENSION(:,:), POINTER :: max_subs
            !!! maximum of extension of subs (single)
            !!! Note: 1st index: x,y,(z), 2nd: subID
            REAL(ppm_kind_double), DIMENSION(:,:), POINTER :: min_subd
            !!! mimimum of extension of subs (double)
            !!! Note: 1st index: x,y,(z), 2nd: subID
            REAL(ppm_kind_double), DIMENSION(:,:), POINTER :: max_subd
            !!! maximum of extension of subs (double)
            !!! Note: 1st index: x,y,(z), 2nd: subID

            REAL(ppm_kind_single), DIMENSION(:  ), POINTER :: sub_costs
            !!! estimated cost associated with subdomains (single). Index: sub-ID.
            REAL(ppm_kind_double), DIMENSION(:  ), POINTER :: sub_costd
            !!! estimated cost associated with subdomains (double). Index: sub-ID.

            INTEGER              , DIMENSION(:  ), POINTER :: sub2proc
            !!! subdomain to processor assignment. index: subID (global)

            INTEGER                                        :: nsublist
            !!! number of subs on the current processor.

            INTEGER              , DIMENSION(:  ), POINTER :: isublist
            !!! list of subs of the current processor. 1st index: local sub
            !!! number.

            INTEGER              , DIMENSION(:,:), POINTER :: subs_bc
            !!! boundary conditions on a sub:
            !!!
            !!! - west  : 1
            !!! - east  : 2
            !!! - south : 3
            !!! - north : 4
            !!! - bottom: 5
            !!! - top   : 6
            !!!
            !!! index 1: the index of the 4 or 6 faces in 2 and 3 D
            !!! index 2: the global sub id
            !!!
            !!! states:
            !!!
            !!! - value: 0 the face is internal
            !!! - value: 1 otherwise

            INTEGER              , DIMENSION(:,:), POINTER :: ineighsubs
            !!! list of neighboring subs of all local subs.
            !!! - index 1: neighbor index
            !!! - index 2: sub id (local index, not global ID!)

            INTEGER              , DIMENSION(:  ), POINTER :: nneighsubs
            !!! number of neighboring subs of all local subs.
            !!!
            !!! index 1: sub id (local index, not global ID!)
            INTEGER              , DIMENSION(:  ), POINTER :: ineighproc
            !!! list of neighboring processors. Index 1: neighbor index
            INTEGER                                        :: nneighproc
            !!! number of neighboring processors.
            LOGICAL                                        :: isoptimized
            !!! has optimal communication sequence already been determined for
            !!! this topology?
            INTEGER                                        :: ncommseq
            !!! number of communication rounds needed for partial mapping
            INTEGER              , DIMENSION(:  ), POINTER :: icommseq
            !!! optimal communication sequence for this processor. 1st index:
            !!! communication round
            INTEGER                                        :: max_meshid
            !!! Number of meshes defined on this topology

            TYPE(ppm_t_equi_mesh), DIMENSION(:  ), POINTER :: mesh
            !!! List of meshes defined on this topology. Index: meshid
         END TYPE ppm_t_topo

         !----------------------------------------------------------------------
         ! Wrapper type to be able to have a pointer array to hold topologies
         !----------------------------------------------------------------------
         TYPE ppm_ptr_t_topo
             TYPE(ppm_t_topo), POINTER  :: t
         END TYPE ppm_ptr_t_topo

         !----------------------------------------------------------------------
         !  Operator interfaces
         !----------------------------------------------------------------------
!         INTERFACE OPERATOR (=)
!             SUBROUTINE ppm_topo_copy(t1,t2,info)
!                 IMPLICIT NONE
!                 TYPE(ppm_t_topo), INTENT(IN   ) :: t1
!                 TYPE(ppm_t_topo), INTENT(  OUT) :: t2
!                 INTEGER         , INTENT(  OUT) :: info
!             END SUBROUTINE
!         END INTERFACE
!         INTERFACE OPERATOR (=)
!             SUBROUTINE ppm_mesh_copy(m1,m2,info)
!                 IMPLICIT NONE
!                 TYPE(ppm_t_equi_mesh), INTENT(IN   ) :: m1
!                 TYPE(ppm_t_equi_mesh), INTENT(  OUT) :: m2
!                 INTEGER              , INTENT(  OUT) :: info
!             END SUBROUTINE
!         END INTERFACE

      END MODULE ppm_module_typedef

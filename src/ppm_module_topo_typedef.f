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

      MODULE ppm_module_topo_typedef
      !!! This module contains the definition of some ppm data types
      !!! and stores them.
      !!! NOTE: they should perhaps be stored in ppm_module_data instead, but
      !!! it is then harder to avoid module circular dependency.

      !----------------------------------------------------------------------
      !  Modules
      !----------------------------------------------------------------------
      USE ppm_module_data

      !----------------------------------------------------------------------
      ! Data types
      !----------------------------------------------------------------------
      !TYPE ppm_t_patch
          !!!! The patch type

          !INTEGER                                       :: topoid
          !!!! id of the topology on which this patch belongs
          !INTEGER              , DIMENSION(:  ),POINTER :: npatch => NULL()
          !!!! number of patches for each subdomain

          !INTEGER              , DIMENSION(:,:),POINTER :: istart => NULL()
          !!!! Starting indices of the mesh of this patch in the mesh for this sub
          !INTEGER              , DIMENSION(:,:),POINTER :: iend => NULL()
          !!!! Ending indices of the mesh of this patch in the mesh for this sub
      !END TYPE ppm_t_patch

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

          REAL(ppm_kind_single), DIMENSION(:), POINTER :: min_physs => NULL()
          !!! minimum of physical extend of the computational domain (single)
          !!! Note: first index is ppm_dim
          REAL(ppm_kind_single), DIMENSION(:), POINTER :: max_physs => NULL()
          !!! maximum of physical extend of the computational domain (single)
          !!! Note: first index is ppm_dim
          REAL(ppm_kind_double), DIMENSION(:), POINTER :: min_physd => NULL()
          !!! minimum of physical extend of the computational domain (double)
          !!! Note: first index is ppm_dim
          REAL(ppm_kind_double), DIMENSION(:), POINTER :: max_physd => NULL()
          !!! maximum of physical extend of the computational domain (double)
          !!! Note: first index is ppm_dim

          INTEGER              , DIMENSION(:  ), POINTER :: bcdef => NULL()
          !!! boundary conditions for the topology
          !!! Note: first index is 1-6 (each of the faces)

          INTEGER                                        :: nsubs
          !!! total number of subs on all processors.

          REAL(ppm_kind_single), DIMENSION(:,:), POINTER :: min_subs => NULL()
          !!! minimum of extension of subs (single)
          !!! Note: 1st index: x,y,(z), 2nd: global subID
          REAL(ppm_kind_single), DIMENSION(:,:), POINTER :: max_subs => NULL()
          !!! maximum of extension of subs (single)
          !!! Note: 1st index: x,y,(z), 2nd: global subID
          REAL(ppm_kind_double), DIMENSION(:,:), POINTER :: min_subd => NULL()
          !!! minimum of extension of subs (double)
          !!! Note: 1st index: x,y,(z), 2nd: global subID
          REAL(ppm_kind_double), DIMENSION(:,:), POINTER :: max_subd => NULL()
          !!! maximum of extension of subs (double)
          !!! Note: 1st index: x,y,(z), 2nd: global subID

          REAL(ppm_kind_single), DIMENSION(:  ),POINTER :: sub_costs => NULL()
          !!! estimated cost associated with subdomains (single). Index: sub-ID.
          REAL(ppm_kind_double), DIMENSION(:  ),POINTER :: sub_costd => NULL()
          !!! estimated cost associated with subdomains (double). Index: sub-ID.

          INTEGER              , DIMENSION(:  ),POINTER :: sub2proc => NULL()
          !!! subdomain to processor assignment. index: subID (global)

          INTEGER                                       :: nsublist
          !!! number of subs on the current processor.

          INTEGER              , DIMENSION(:  ),POINTER :: isublist => NULL()
          !!! list of global sub IDs of the current processor.
          !!! 1st index: local sub index.
          !!! Given the local sub ID, you get the global ID

          INTEGER              , DIMENSION(:,:),POINTER :: subs_bc => NULL()
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
          INTEGER            , DIMENSION(:,:), POINTER :: ineighsubs => NULL()
          !!! list of neighboring subs (global IDs) of all local subs.
          !!! - index 1: neighbor index
          !!! - index 2: sub id (local index, not global ID!)

          INTEGER            , DIMENSION(:  ), POINTER :: nneighsubs => NULL()
          !!! number of neighboring subs of all local subs.
          !!!
          !!! index 1: sub id (local index, not global ID!)
          INTEGER            , DIMENSION(:  ), POINTER :: ineighproc => NULL()
          !!! list of neighboring processors. Index 1: neighbor index
          INTEGER                                      :: nneighproc
          !!! number of neighboring processors.
          LOGICAL                                      :: isoptimized
          !!! has optimal communication sequence already been determined for
          !!! this topology?
          INTEGER                                      :: ncommseq
          !!! number of communication rounds needed for partial mapping
          INTEGER            , DIMENSION(:  ), POINTER :: icommseq => NULL()
          !!! optimal communication sequence for this processor. 1st index:
          !!! communication round
          !          INTEGER                                        :: max_meshid
          !          !!! Number of meshes defined on this topology

          !TYPE(ppm_t_patch), DIMENSION(:),POINTER :: patches => NULL()
          !!!! List of patches data structures on this topology.

          REAL(ppm_kind_single)                        :: ghostsizes
          !!! max ghostsize width used when creating this topology (single)
          REAL(ppm_kind_double)                        :: ghostsized
          !!! max ghostsize width used when creating this topology (double)
          !!! using a larger cutoff when calling routines such as ghost_get
          !!! should ideally raise a warning
      END TYPE ppm_t_topo

      !----------------------------------------------------------------------
      ! Wrapper type to be able to have a pointer array to hold topologies
      !----------------------------------------------------------------------
      TYPE ppm_ptr_t_topo
          TYPE(ppm_t_topo), POINTER  :: t => NULL()
      END TYPE ppm_ptr_t_topo


      !----------------------------------------------------------------------
      ! Pointer to cell list (needed to make lists of cell lists)
      !----------------------------------------------------------------------
      TYPE ppm_t_clist
          !!! Cell list data structure
          INTEGER, DIMENSION(:), POINTER    :: nm  => NULL()
          !!! Number of cells in x,y,(z) direction (including the ghosts 
          !!! cells) in each subdomain. 
          INTEGER, DIMENSION(:), POINTER    :: lpdx => NULL()
          !!! particle index list
          INTEGER, DIMENSION(:), POINTER    :: lhbx => NULL()
          !!! first particle in each cell
      END TYPE


      !----------------------------------------------------------------------
      ! DATA STORAGE for the topologies
      !----------------------------------------------------------------------
      TYPE(ppm_ptr_t_topo), DIMENSION(:), POINTER :: ppm_topo => null()
      !!! the ppm topologies array

      INTEGER :: ppm_next_avail_topo
      !!! id of the next available topology to be used by
      !!! ppm_topo_alloc.
      !!! 
      !!! at initialization this is set to ppm_param_undefined              +
      !!! if it points within (1,size(ppm_topo)) this slot is (re)used to
      !!! store the new topology, if it is > size(ppm_topo) then the ppm_topo
      !!! array must be extended



      END MODULE ppm_module_topo_typedef

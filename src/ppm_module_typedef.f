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
      !!! This module contains the definition of some ppm data types
      !!! and stores them.
      !!! NOTE: they should perhaps be stored in ppm_module_data instead, but
      !!! it is then harder to avoid module circular dependency.

         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_mesh_typedef
         USE ppm_module_particles_typedef
         USE ppm_module_topo_typedef

         !----------------------------------------------------------------------
         ! Data types
         !----------------------------------------------------------------------
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
         ! Global variables and parameters
         !----------------------------------------------------------------------

         !----------------------------------------------------------------------
         ! Topologies
         !----------------------------------------------------------------------
         TYPE(ppm_ptr_t_topo_s), DIMENSION(:), POINTER :: ppm_topo_s => NULL()
         !!! the PPM topologies array (single precision)
         TYPE(ppm_ptr_t_topo_d), DIMENSION(:), POINTER :: ppm_topo_d => NULL()
         !!! the PPM topologies array (double precision)

         INTEGER :: ppm_next_avail_topo
         !!! ID of the next available topology to be used by
         !!! ppm_topo_alloc.
         !!!
         !!! At initialization this is set to ppm_param_undefined              +
         !!! If it points within (1,SIZE(ppm_topo)) this slot is (re)used to
         !!! store the new topology, if it is > SIZE(ppm_topo) then the ppm_topo
         !!! array must be extended
      END MODULE ppm_module_typedef

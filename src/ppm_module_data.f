      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_data
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
         !  mapping
         !----------------------------------------------------------------------
         INTEGER                        :: ppm_map_type 
         INTEGER                        :: ppm_nsendlist
         INTEGER                        :: ppm_nrecvlist
         INTEGER, DIMENSION(:), POINTER :: ppm_isendlist => NULL()
         INTEGER, DIMENSION(:), POINTER :: ppm_irecvlist => NULL()

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
         TYPE(ppm_ptr_t_topo), DIMENSION(:), POINTER :: ppm_topo => NULL()
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
         REAL(ppm_kind_double),DIMENSION(:),POINTER :: ppm_proc_speed => NULL()

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

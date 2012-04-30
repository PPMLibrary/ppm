      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_data_loadbal
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
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
     
      MODULE ppm_module_data_loadbal
      !!! This module holds data used by the load balancing routines.
      !!!
      !!! [NOTE]
      !!! The variables declared in this modules should not be directly accessed
      !!! by the user. They are managed interally by the library.
         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_data, ONLY: ppm_kind_double
         PRIVATE :: ppm_kind_double

         !----------------------------------------------------------------------
         !  Timing and load statistics
         !----------------------------------------------------------------------
         ! estimated cost of redefining the topology
         REAL(ppm_kind_double) :: ppm_loadbal_decomp_cost = 0.0_ppm_kind_double
         ! counter of how many times decomp cost was measured (N)
         INTEGER               :: ppm_loadbal_dcn = 0
         ! running sum of (maxtime - avgtime)
         REAL(ppm_kind_double) :: ppm_loadbal_runsum = 0.0_ppm_kind_double
         ! old value of the SAR function
         REAL(ppm_kind_double) :: ppm_loadbal_old_sar = -1.0_ppm_kind_double

      END MODULE ppm_module_data_loadbal

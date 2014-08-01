      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                ppm_module_topo
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

      MODULE ppm_module_topo
      !!! This module contains all user-callable routines
      !!! needed to create and manage PPM topologies.
         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_topo_typedef, ONLY : ppm_t_topo,ppm_topo
         USE ppm_module_mktopo,       ONLY : ppm_mktopo
         USE ppm_module_topo_check,   ONLY : ppm_topo_check
         USE ppm_module_mesh_define,  ONLY : ppm_mesh_define
         USE ppm_module_scale_domain, ONLY : ppm_scale_domain
         USE ppm_module_topo_get,     ONLY : ppm_topo_get, &
         &   ppm_topo_get_decomp

      END MODULE ppm_module_topo

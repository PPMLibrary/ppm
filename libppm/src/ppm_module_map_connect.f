      !--*- f90 -*--------------------------------------------------------------
      !  Module       :            ppm_module_map_connect_distrib
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
     
      MODULE ppm_module_map_connect
      !!! This module provides the mapping routines for connections
         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:)  , POINTER :: id_send,id_temp,id_inv,csend,crecv
         INTEGER, DIMENSION(:)  , POINTER :: sendbuffer,recvbuffer,tempbuffer
         INTEGER, DIMENSION(:,:), POINTER :: cd_local,psend

         PRIVATE :: id_send,id_temp,id_inv,csend,crecv
         PRIVATE :: sendbuffer,recvbuffer,tempbuffer
         PRIVATE :: psend,cd_local

         !----------------------------------------------------------------------
         !  Define inferface to ppm_map_connect_distrib
         !----------------------------------------------------------------------
         INTERFACE ppm_map_connect_distrib
            MODULE PROCEDURE ppm_map_connect_distrib
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to ppm_map_connect_prune
         !----------------------------------------------------------------------
         INTERFACE ppm_map_connect_prune
            MODULE PROCEDURE ppm_map_connect_prune
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define inteface to ppm_map_connect_send
         !----------------------------------------------------------------------
         INTERFACE ppm_map_connect_send
            MODULE PROCEDURE ppm_map_connect_send
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "map/ppm_map_connect_distrib.f"

#include "map/ppm_map_connect_prune.f"

#include "map/ppm_map_connect_send.f"

      END MODULE ppm_module_map_connect

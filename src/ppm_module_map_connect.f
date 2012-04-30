      !--*- f90 -*--------------------------------------------------------------
      !  Module       :            ppm_module_map_connect_distrib
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
     
      MODULE ppm_module_map_connect
      !!! This module provides the mapping routines for connections
         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:)  , POINTER :: id_send => NULL()
         INTEGER, DIMENSION(:)  , POINTER :: id_temp => NULL()
         INTEGER, DIMENSION(:)  , POINTER :: id_inv  => NULL()
         INTEGER, DIMENSION(:)  , POINTER :: csend   => NULL()
         INTEGER, DIMENSION(:)  , POINTER :: crecv   => NULL()
         INTEGER, DIMENSION(:)  , POINTER :: sendbuffer => NULL()
         INTEGER, DIMENSION(:)  , POINTER :: recvbuffer => NULL()
         INTEGER, DIMENSION(:)  , POINTER :: tempbuffer => NULL()
         INTEGER, DIMENSION(:,:), POINTER :: cd_local   => NULL()
         INTEGER, DIMENSION(:,:), POINTER :: psend      => NULL()

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

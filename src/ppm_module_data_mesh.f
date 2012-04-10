      !--*- f90 -*--------------------------------------------------------------
      !  Module       :               ppm_module_data_mesh
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

      MODULE ppm_module_data_mesh
      !!! This module holds mesh mapping, receive buffers and
      !!! mesh ghost layer mapping buffers.
      !!!
      !!! [NOTE]
      !!! The variables declared in this module should not be accessed by the
      !!! PPM client developer. They are managed interally by the library.

         !----------------------------------------------------------------------
         !  Mesh mapping, send and receive lists
         !----------------------------------------------------------------------

         INTEGER, DIMENSION(:  ), POINTER :: ppm_mesh_isendfromsub => NULL()
         !!! list of source subs to send from local processor (local sub number
         !!! on source processor)
         INTEGER, DIMENSION(:,:), POINTER :: ppm_mesh_isendblkstart => NULL()
         ! start (lower-left corner) of mesh block to be sent in GLOBAL
         ! mesh coordinates. First index: x,y[,z], 2nd: isendlist
         INTEGER, DIMENSION(:,:), POINTER :: ppm_mesh_isendpatchid => NULL()
         !!! list of source patch ids to send from local processor 
         INTEGER, DIMENSION(:,:), POINTER :: ppm_mesh_isendblksize => NULL()
         ! size (in grid points) of blocks to be sent
         INTEGER, DIMENSION(:  ), POINTER :: ppm_mesh_irecvtosub => NULL()
         ! list of destination subs to recv to on local processors (local sub
         ! number on destination processor)
         INTEGER, DIMENSION(:,:), POINTER :: ppm_mesh_irecvpatchid => NULL()
         !!! list of source patch ids to receive from local processor 
         INTEGER, DIMENSION(:,:), POINTER :: ppm_mesh_irecvblkstart => NULL()
         ! start (lower-left corner) of mesh block to be recvd in GLOBAL
         ! mesh coordinates. First index: x,y[,z], 2nd: isendlist
         INTEGER, DIMENSION(:,:), POINTER :: ppm_mesh_irecvblksize => NULL()
         ! size (in grid points) of blocks to be recvd



      END MODULE ppm_module_data_mesh

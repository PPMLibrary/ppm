      !--*- f90 -*--------------------------------------------------------------
      !  Module       :           ppm_module_map_field_ghost
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_map_field_ghost
      !!! This module contains interfaces to the field ghost mapping routines
      !!! and all data structures and definitions that
      !!! are `PRIVATE` to the mesh routines.
      !!!
      !!! [NOTE]
      !!! The terminology distinguishes between meshes and fields
      !!! (the data living on the meshes). Several fields can use the
      !!! same mesh. Meshes are defined as per-topology ppm-internal types,
      !!! whereas fields are user-provided arrays.
         !----------------------------------------------------------------------
         !  Work memory
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:  ), POINTER :: isendfromsub,isendtosub
         INTEGER, DIMENSION(:  ), POINTER :: sendbuf,recvbuf
         INTEGER, DIMENSION(:,:), POINTER :: isendblkstart,isendblksize,ioffset
         ! sorted (according to proc-proc interaction order) offset list)
         INTEGER, DIMENSION(:,:), POINTER :: mesh_ghost_offset

         PRIVATE :: isendfromsub,isendtosub,sendbuf,recvbuf,isendblkstart
         PRIVATE :: isendblksize,ioffset,mesh_ghost_offset

         !----------------------------------------------------------------------
         !  Define interface to ppm_map_field_ghost_init
         !----------------------------------------------------------------------
         INTERFACE ppm_map_field_ghost_init
             MODULE PROCEDURE ppm_map_field_ghost_init
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to ppm_map_field_ghost_get
         !----------------------------------------------------------------------
         INTERFACE ppm_map_field_ghost_get
             MODULE PROCEDURE ppm_map_field_ghost_get
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to ppm_map_field_ghost_put
         !----------------------------------------------------------------------
         INTERFACE ppm_map_field_ghost_put
             MODULE PROCEDURE ppm_map_field_ghost_put
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "map/ppm_map_field_ghost_init.f"

#include "map/ppm_map_field_ghost_get.f"

#include "map/ppm_map_field_ghost_put.f"

      END MODULE ppm_module_map_field_ghost

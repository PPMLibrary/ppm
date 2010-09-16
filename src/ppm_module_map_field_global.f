      !--*- f90 -*--------------------------------------------------------------
      !  Module       :            ppm_module_map_field_global
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_map_field_global
      !!! This module contains the interface and the needed work arrays for the
      !!! field global mappings.
      !!!
      !!! [NOTE]
      !!! The terminology distinguishes between meshes and fields
      !!! (the data living on the meshes). Several fields can use the
      !!! same mesh. Meshes are defined as ppm-internal TYPES, whereas
      !!! fields are user-provided arrays.
         !----------------------------------------------------------------------
         !  Work memory
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:  ), POINTER :: isendfromsub,isendtosub
         INTEGER, DIMENSION(:,:), POINTER :: isendblkstart,isendblksize,ioffset
         INTEGER, DIMENSION(:  ), POINTER :: irecvfromsub,irecvtosub
         INTEGER, DIMENSION(:,:), POINTER :: irecvblkstart,irecvblksize

         PRIVATE :: isendfromsub,isendtosub,isendblkstart,isendblksize
         PRIVATE :: ioffset,irecvfromsub,irecvtosub,irecvblkstart,irecvblksize

         !----------------------------------------------------------------------
         !  Define interface to ppm_map_field_global
         !----------------------------------------------------------------------
         INTERFACE ppm_map_field_global
             MODULE PROCEDURE ppm_map_field_global
         END INTERFACE

         INTERFACE ppm_map_field_global_symm
            MODULE PROCEDURE ppm_map_field_global_symm
         END INTERFACE

         INTERFACE ppm_map_field_global_useperiod
            MODULE PROCEDURE ppm_map_field_global_useperiod
         END INTERFACE

         INTERFACE ppm_map_field_global_useperiod_store
            MODULE PROCEDURE ppm_map_field_global_useperiod_store
         END INTERFACE

         INTERFACE ppm_map_field_globalstored
            MODULE PROCEDURE ppm_map_field_globalstored
         END INTERFACE
         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "map/ppm_map_field_global.f"

#include "map/ppm_map_field_global_symm.f"

#include "map/ppm_map_field_global_useperiod.f"

#include "map/ppm_map_field_global_useperiod_store.f"

#include "map/ppm_map_field_globalstored.f"

      END MODULE ppm_module_map_field_global

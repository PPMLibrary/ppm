      !--*- f90 -*--------------------------------------------------------------
      !  Module       :               ppm_module_mesh_define
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2

      MODULE ppm_module_mesh_define
      !!! This module contains the Interface to `ppm_mesh_define`
      !!!
      !!! [NOTE]
      !!! The terminology distinguishes between meshes and fields
      !!! (the data living on the meshes). Several fields can use the
      !!! same mesh. Meshes are defined as ppm-internal TYPES, whereas
      !!! fields are user-provided arrays.
         !----------------------------------------------------------------------
         !  Define interfaces to mesh creation routine
         !----------------------------------------------------------------------
         INTERFACE ppm_mesh_define
            MODULE PROCEDURE ppm_mesh_define
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "mesh/ppm_mesh_define.f"


      END MODULE ppm_module_mesh_define

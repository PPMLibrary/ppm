      !--*- f90 -*--------------------------------------------------------------
      !  Module       :               ppm_module_mesh_derive
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_mesh_derive
      !!! This module contains the interface to `ppm_mesh_derive`.
      !!!
      !!! [NOTE]
      !!! The terminology distinguishes between meshes and fields
      !!! (the data living on the meshes). Several fields can use the
      !!! same mesh. Meshes are defined as ppm-internal TYPES, whereas
      !!! fields are user-provided arrays.
         !----------------------------------------------------------------------
         !  Define interface to ppm_mesh_derive
         !----------------------------------------------------------------------
         INTERFACE ppm_mesh_derive
             MODULE PROCEDURE ppm_mesh_derive
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "mesh/ppm_mesh_derive.f"

      END MODULE ppm_module_mesh_derive

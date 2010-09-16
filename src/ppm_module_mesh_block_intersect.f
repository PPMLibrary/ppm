      !--*- f90 -*--------------------------------------------------------------
      !  Module       :          ppm_module_mesh_block_intersect
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_mesh_block_intersect
      !!! This module contains the Interface to ppm_mesh_block_intersect
      !!!
      !!! [NOTE]
      !!! The terminology distinguishes between meshes and fields
      !!! (the data living on the meshes). Several fields can use the
      !!! same mesh. Meshes are defined as ppm-internal TYPES, whereas
      !!! fields are user-provided arrays.
         !----------------------------------------------------------------------
         !  Define interface to ppm_mesh_block_intersect
         !----------------------------------------------------------------------
         INTERFACE ppm_mesh_block_intersect
             MODULE PROCEDURE ppm_mesh_block_intersect
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "mesh/ppm_mesh_block_intersect.f"

      END MODULE ppm_module_mesh_block_intersect

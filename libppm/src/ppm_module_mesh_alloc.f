      !--*- f90 -*--------------------------------------------------------------
      !  Module       :             ppm_module_mesh_alloc
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_mesh_alloc
      !!! This module contains the Interface to ppm_mesh_alloc.
      !!!
      !!! [NOTE]
      !!! The terminology distinguishes between meshes and fields
      !!! (the data living on the meshes). Several fields can use the
      !!! same mesh. Meshes are defined as ppm-internal TYPES, whereas
      !!! fields are user-provided arrays.
         !----------------------------------------------------------------------
         !  Define interface to ppm_mesh_alloc
         !----------------------------------------------------------------------
         INTERFACE ppm_mesh_alloc
             MODULE PROCEDURE ppm_mesh_alloc_equi
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "mesh/ppm_mesh_alloc_equi.f"

      END MODULE ppm_module_mesh_alloc

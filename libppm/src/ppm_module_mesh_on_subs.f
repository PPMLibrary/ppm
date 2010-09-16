      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_mesh_on_subs
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

      MODULE ppm_module_mesh_on_subs
      !!! This module contains the interface to `ppm_mesh_on_subs`.
      !!!
      !!! [NOTE]
      !!! The terminology distinguishes between meshes and fields
      !!! (the data living on the meshes). Several fields can use the
      !!! same mesh. Meshes are defined as ppm-internal TYPES, whereas
      !!! fields are user-provided arrays.
         !----------------------------------------------------------------------
         !  Define interface to ppm_mesh_on_subs
         !----------------------------------------------------------------------
         INTERFACE ppm_mesh_on_subs
             MODULE PROCEDURE ppm_mesh_on_subs_s
             MODULE PROCEDURE ppm_mesh_on_subs_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "mesh/ppm_mesh_on_subs.f"
#undef __KIND
         
#define __KIND __DOUBLE_PRECISION
#include "mesh/ppm_mesh_on_subs.f"
#undef __KIND

      END MODULE ppm_module_mesh_on_subs

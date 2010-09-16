      !--*- f90 -*--------------------------------------------------------------
      !  Module       :             ppm_module_mesh_finalize
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
#define __INTEGER                  3
#define __LOGICAL                  4
#define __SINGLE_PRECISION_COMPLEX 5
#define __DOUBLE_PRECISION_COMPLEX 6
#define __2D                       7
#define __3D                       8
#define __SFIELD                   9
#define __VFIELD                   10

      MODULE ppm_module_mesh_finalize
      !!! This module contains the interface and the different implementations
      !!! of ppm_mesh_finalize.
      !!!
      !!! [NOTE]
      !!! The terminology distinguishes between meshes and fields
      !!! (the data living on the meshes). Several fields can use the
      !!! same mesh. Meshes are defined as ppm-internal TYPES, whereas
      !!! fields are user-provided arrays.
         !----------------------------------------------------------------------
         !  Define interface to ppm_mesh_finalize
         !----------------------------------------------------------------------
         INTERFACE ppm_mesh_finalize
             MODULE PROCEDURE ppm_mesh_finalize
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "mesh/ppm_mesh_finalize.f"

      END MODULE ppm_module_mesh_finalize

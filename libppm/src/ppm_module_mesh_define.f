      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_mesh_define
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all data structures and
      !                 definitions that are PRIVATE to the mesh routines.
      !                 It also included those routines and provides
      !                 INTERFACEs.
      !                
      !  Remarks      : The terminology distinguishes between meshes and
      !                 fields (the data living on the meshes). Several
      !                 fields can use the same mesh. Meshes are defined as
      !                 ppm-internal TYPES, whereas fields are
      !                 user-provided arrays.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_mesh_define.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:29:58  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2

      MODULE ppm_module_mesh_define

         !----------------------------------------------------------------------
         !  Define interfaces to mesh creation routine
         !----------------------------------------------------------------------
         INTERFACE ppm_mesh_define
            MODULE PROCEDURE ppm_mesh_define_s
            MODULE PROCEDURE ppm_mesh_define_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_mesh_define.f"
#undef __KIND
         
#define __KIND __DOUBLE_PRECISION
#include "ppm_mesh_define.f"
#undef __KIND

      END MODULE ppm_module_mesh_define

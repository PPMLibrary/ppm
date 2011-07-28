      !-------------------------------------------------------------------------
      !  Module       :           ppm_module_map_field_ghost_put
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
      !  $Log: ppm_module_map_field_ghost_put.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:59  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:29:45  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_map_field_ghost_put

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

#include "ppm_map_field_ghost_put.f"

      END MODULE ppm_module_map_field_ghost_put

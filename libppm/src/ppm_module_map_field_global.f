      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_map_field_global
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
      !  $Log: ppm_module_map_field_global.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:59  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/11/11 15:26:18  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.1  2004/07/26 07:29:46  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_map_field_global

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
             ! We do not need different accuracy versions since field
             ! indices are always INTEGER. Moreover, there is no partial
             ! mapping since mesh points cannot move. Every map is thus
             ! going onto a new topology.
             MODULE PROCEDURE ppm_map_field_global
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_map_field_global.f"

      END MODULE ppm_module_map_field_global

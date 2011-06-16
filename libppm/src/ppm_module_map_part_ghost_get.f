      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_map_part_ghost_get
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the mapping
      !                 routines. 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_map_part_ghost_get.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/11/11 15:26:19  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.1  2004/07/26 07:29:52  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2

      MODULE ppm_module_map_part_ghost_get

         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:), POINTER :: ilist1,ilist2,ighost
         LOGICAL, DIMENSION(:), POINTER :: lghost

         PRIVATE :: ilist1,ilist2,ighost,lghost

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_ghost_get
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_ghost_get
            MODULE PROCEDURE ppm_map_part_ghost_get_d
            MODULE PROCEDURE ppm_map_part_ghost_get_s
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_map_part_ghost_get.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_map_part_ghost_get.f"
#undef __KIND

      END MODULE ppm_module_map_part_ghost_get

      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_map_part_ring
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
      !  $Log: ppm_module_map_part_ring.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/11/11 15:26:19  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.1  2004/07/26 07:29:55  ivos
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

      MODULE ppm_module_map_part_ring

         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:), POINTER :: plist_des,plist_act,plist_exc
         INTEGER, DIMENSION(:), POINTER :: srlist1,srlist2

         PRIVATE :: plist_des,plist_act,plist_exc,srlist1,srlist2

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_ring
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_ring
            MODULE PROCEDURE ppm_map_part_ring_d
            MODULE PROCEDURE ppm_map_part_ring_s
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_map_part_ring.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_map_part_ring.f"
#undef __KIND

      END MODULE ppm_module_map_part_ring

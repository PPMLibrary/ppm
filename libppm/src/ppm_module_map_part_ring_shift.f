      !-------------------------------------------------------------------------
      !  Module       :           ppm_module_map_part_ring_shift
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
      !  $Log: ppm_module_map_part_ring_shift.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:29:56  ivos
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

      MODULE ppm_module_map_part_ring_shift

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_ring_shift
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_ring_shift
            MODULE PROCEDURE ppm_map_part_ring_shift_s
            MODULE PROCEDURE ppm_map_part_ring_shift_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_map_part_ring_shift.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_map_part_ring_shift.f"
#undef __KIND

      END MODULE ppm_module_map_part_ring_shift

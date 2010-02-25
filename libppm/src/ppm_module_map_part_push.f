      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_map_part_push
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
      !  $Log: ppm_module_map_part_push.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:29:54  ivos
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
#define __INTEGER                  3
#define __LOGICAL                  4
#define __SINGLE_PRECISION_COMPLEX 5
#define __DOUBLE_PRECISION_COMPLEX 6

      MODULE ppm_module_map_part_push

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_push
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_push
            ! scalar (1d) particle data
            MODULE PROCEDURE ppm_map_part_push_1dd
            MODULE PROCEDURE ppm_map_part_push_1ds
            MODULE PROCEDURE ppm_map_part_push_1di
            MODULE PROCEDURE ppm_map_part_push_1dl
            MODULE PROCEDURE ppm_map_part_push_1ddc
            MODULE PROCEDURE ppm_map_part_push_1dsc

            ! vector (2d) particle data
            MODULE PROCEDURE ppm_map_part_push_2dd
            MODULE PROCEDURE ppm_map_part_push_2ds
            MODULE PROCEDURE ppm_map_part_push_2di
            MODULE PROCEDURE ppm_map_part_push_2dl
            MODULE PROCEDURE ppm_map_part_push_2ddc
            MODULE PROCEDURE ppm_map_part_push_2dsc
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_map_part_push.f"
#undef __KIND
#define __KIND __INTEGER
#include "ppm_map_part_push.f"
#undef __KIND
#define __KIND __LOGICAL
#include "ppm_map_part_push.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_map_part_push.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_map_part_push.f"
#undef __KIND
#define __KIND __INTEGER
#include "ppm_map_part_push.f"
#undef __KIND
#define __KIND __LOGICAL
#include "ppm_map_part_push.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_map_part_push.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_map_part_push.f"
#undef __KIND
#undef __DIM

      END MODULE ppm_module_map_part_push

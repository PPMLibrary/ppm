      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_map_part_ghost
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
      !  $Log: ppm_module_map_part_ghost.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 07:29:51  ivos
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

      MODULE ppm_module_map_part_ghost

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_ghost
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_ghost
            ! scalar (1d) particle data
            MODULE PROCEDURE ppm_map_part_ghost_1dd
            MODULE PROCEDURE ppm_map_part_ghost_1ds
            MODULE PROCEDURE ppm_map_part_ghost_1ddc
            MODULE PROCEDURE ppm_map_part_ghost_1dsc
            MODULE PROCEDURE ppm_map_part_ghost_1dis
            MODULE PROCEDURE ppm_map_part_ghost_1did
            MODULE PROCEDURE ppm_map_part_ghost_1dls
            MODULE PROCEDURE ppm_map_part_ghost_1dld

            ! vector (2d) particle data
            MODULE PROCEDURE ppm_map_part_ghost_2dd
            MODULE PROCEDURE ppm_map_part_ghost_2ds
            MODULE PROCEDURE ppm_map_part_ghost_2ddc
            MODULE PROCEDURE ppm_map_part_ghost_2dsc
            MODULE PROCEDURE ppm_map_part_ghost_2dis
            MODULE PROCEDURE ppm_map_part_ghost_2did
            MODULE PROCEDURE ppm_map_part_ghost_2dls
            MODULE PROCEDURE ppm_map_part_ghost_2dld
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "ppm_map_part_ghost.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_map_part_ghost.f"
#undef __KIND
#define __KIND __INTEGER
#define __KIND_AUX __DOUBLE_PRECISION
#include "ppm_map_part_ghost.f"
#undef __KIND_AUX
#define __KIND_AUX __SINGLE_PRECISION
#include "ppm_map_part_ghost.f"
#undef __KIND_AUX
#undef __KIND
#define __KIND __LOGICAL
#define __KIND_AUX __DOUBLE_PRECISION
#include "ppm_map_part_ghost.f"
#undef __KIND_AUX
#define __KIND_AUX __SINGLE_PRECISION
#include "ppm_map_part_ghost.f"
#undef __KIND_AUX
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_map_part_ghost.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_map_part_ghost.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "ppm_map_part_ghost.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_map_part_ghost.f"
#undef __KIND
#define __KIND __INTEGER
#define __KIND_AUX __DOUBLE_PRECISION
#include "ppm_map_part_ghost.f"
#undef __KIND_AUX
#define __KIND_AUX __SINGLE_PRECISION
#include "ppm_map_part_ghost.f"
#undef __KIND_AUX
#undef __KIND
#define __KIND __LOGICAL
#define __KIND_AUX __DOUBLE_PRECISION
#include "ppm_map_part_ghost.f"
#undef __KIND_AUX
#define __KIND_AUX __SINGLE_PRECISION
#include "ppm_map_part_ghost.f"
#undef __KIND_AUX
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_map_part_ghost.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_map_part_ghost.f"
#undef __KIND
#undef __DIM

      END MODULE ppm_module_map_part_ghost

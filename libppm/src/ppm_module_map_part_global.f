      !--*- f90 -*--------------------------------------------------------------
      !  Module       :             ppm_module_map_part_global
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

      MODULE ppm_module_map_part_global
      !!! This module provides the mapping routines for global particle
      !!! mapping
         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:), POINTER :: ilist1,ilist2,part2proc

         PRIVATE :: ilist1,ilist2,part2proc

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_global
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_global
            MODULE PROCEDURE ppm_map_part_global_d
            MODULE PROCEDURE ppm_map_part_global_s
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_remap
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_remap
            MODULE PROCEDURE ppm_map_part_remap_d
            MODULE PROCEDURE ppm_map_part_remap_s
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_global.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_global.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_remap.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_remap.f"
#undef __KIND

      END MODULE ppm_module_map_part_global

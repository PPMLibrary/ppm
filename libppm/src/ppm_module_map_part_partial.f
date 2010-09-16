      !--*- f90 -*--------------------------------------------------------------
      !  Module       :             ppm_module_map_part_partial
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

      MODULE ppm_module_map_part_partial
      !!! This module provides the mapping routines for partial particle
      !!! mapping
         !----------------------------------------------------------------------
         !  Work lists
         !----------------------------------------------------------------------
         INTEGER, DIMENSION(:), POINTER :: ilist1,ilist2,part2proc,ineighsubs

         PRIVATE :: ilist1,ilist2,part2proc,ineighsubs

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_partial
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_partial
            MODULE PROCEDURE ppm_map_part_partial_d
            MODULE PROCEDURE ppm_map_part_partial_s
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "map/ppm_map_part_partial.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "map/ppm_map_part_partial.f"
#undef __KIND

      END MODULE ppm_module_map_part_partial

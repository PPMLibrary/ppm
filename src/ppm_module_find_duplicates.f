      !--*- f90 -*--------------------------------------------------------------
      !  Module       :            ppm_module_find_duplicates
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_find_duplicates
      !!! This module provides the routines to find doublicate entries in 2D
      !!! arrays - callable from the outside.
         !----------------------------------------------------------------------
         !  Define interface to duplicate finder
         !----------------------------------------------------------------------
         INTERFACE ppm_find_duplicates
            MODULE PROCEDURE ppm_find_duplicates_s
            MODULE PROCEDURE ppm_find_duplicates_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "util/ppm_find_duplicates.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "util/ppm_find_duplicates.f"
#undef __KIND

      END MODULE ppm_module_find_duplicates

      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                ppm_module_util_sort
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

      MODULE ppm_module_util_sort
      !!! This module provides the utility sorting routines.
         !----------------------------------------------------------------------
         !  Define interfaces to the sorting routine
         !----------------------------------------------------------------------
         INTERFACE ppm_util_sort3d
            MODULE PROCEDURE ppm_util_sort3d_s
            MODULE PROCEDURE ppm_util_sort3d_d
         END INTERFACE
         INTERFACE ppm_util_sort2d
            MODULE PROCEDURE ppm_util_sort2d_s
            MODULE PROCEDURE ppm_util_sort2d_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "util/ppm_util_sort3d.f"
#include "util/ppm_util_sort2d.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "util/ppm_util_sort3d.f"
#include "util/ppm_util_sort2d.f"
#undef __KIND

      END MODULE ppm_module_util_sort

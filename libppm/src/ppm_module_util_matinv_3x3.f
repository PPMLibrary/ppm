      !--*- f90 -*--------------------------------------------------------------
      !  Module       :             ppm_module_util_matinv_3x3
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

      MODULE ppm_module_util_matinv_3x3
      !!! This module provides the routines
      !!! that compute inverse of 3x3 matrices.
         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_util_matinv_3x3
            MODULE PROCEDURE ppm_util_matinv_3x3_s
            MODULE PROCEDURE ppm_util_matinv_3x3_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "util/ppm_util_matinv_3x3.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "util/ppm_util_matinv_3x3.f"
#undef __KIND

      END MODULE ppm_module_util_matinv_3x3

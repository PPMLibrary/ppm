      !--*- f90 -*--------------------------------------------------------------
      !  Module       :             ppm_module_util_matinv_2x2
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

      MODULE ppm_module_util_matinv_2x2
      !!! This module provides the routines
      !!! that compute inverse of 2x2 matrices.
         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_util_matinv_2x2
            MODULE PROCEDURE ppm_util_matinv_2x2_s
            MODULE PROCEDURE ppm_util_matinv_2x2_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "util/ppm_util_matinv_2x2.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "util/ppm_util_matinv_2x2.f"
#undef __KIND

      END MODULE ppm_module_util_matinv_2x2

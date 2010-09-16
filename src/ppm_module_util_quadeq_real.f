      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_util_quadeq_real
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

      MODULE ppm_module_util_quadeq_real
      !!! This module provides the routines
      !!! that solve quadratic equations with real roots.
         !----------------------------------------------------------------------
         !  Define interfaces to the main topology routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_util_quadeq_real
            MODULE PROCEDURE ppm_util_quadeq_real_s
            MODULE PROCEDURE ppm_util_quadeq_real_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "util/ppm_util_quadeq_real.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "util/ppm_util_quadeq_real.f"
#undef __KIND

      END MODULE ppm_module_util_quadeq_real

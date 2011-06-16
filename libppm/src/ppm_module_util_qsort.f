      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_util_qsort
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#define __SINGLE_PRECISION           1
#define __DOUBLE_PRECISION           2
#define __INTEGER                    3

      MODULE ppm_module_util_qsort
      !!! This module provides the utility quicksort routines.
         !----------------------------------------------------------------------
         !  Define interface to the list sorting routine
         !----------------------------------------------------------------------
         INTERFACE ppm_util_qsort
            MODULE PROCEDURE ppm_util_qsort_s
            MODULE PROCEDURE ppm_util_qsort_d
            MODULE PROCEDURE ppm_util_qsort_i
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS
#define __KIND __SINGLE_PRECISION
#include "util/ppm_util_qsort.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "util/ppm_util_qsort.f"
#undef __KIND

#define __KIND __INTEGER
#include "util/ppm_util_qsort.f"
#undef __KIND

      END MODULE ppm_module_util_qsort

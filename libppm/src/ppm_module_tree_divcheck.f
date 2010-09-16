      !--*- f90 -*--------------------------------------------------------------
      !  Module       :               ppm_module_tree_divcheck
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

      MODULE ppm_module_tree_divcheck
      !!! This module provides the
      !!! routine that checks if a box is subdivisible.
         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_tree_divcheck
            MODULE PROCEDURE ppm_tree_divcheck_s
            MODULE PROCEDURE ppm_tree_divcheck_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "tree/ppm_tree_divcheck.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "tree/ppm_tree_divcheck.f"
#undef __KIND

      END MODULE ppm_module_tree_divcheck

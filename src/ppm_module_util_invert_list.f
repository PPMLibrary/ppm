      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_util_invert_list
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_util_invert_list
      !!! This module provides the utility routines for inverting lists
         !----------------------------------------------------------------------
         !  Define interface to the list inversion routine
         !----------------------------------------------------------------------
         INTERFACE ppm_util_invert_list
            MODULE PROCEDURE ppm_util_invert_list
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "util/ppm_util_invert_list.f"

      END MODULE ppm_module_util_invert_list
